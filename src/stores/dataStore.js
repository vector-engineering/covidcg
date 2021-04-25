import { observable, action, toJS } from 'mobx';
import { hostname } from '../config';
import {
  processSelectedSnvs,
  processCooccurrenceData,
} from '../utils/snpDataWorkerWrapper';
import { downloadBlobURL } from '../utils/download';
import { intToISO } from '../utils/date';
import { getLocationIdsByNode } from '../utils/location';
import { formatSnv } from '../utils/snpUtils';
import { asyncDataStoreInstance } from '../components/App';
import { rootStoreInstance } from './rootStore';
import {
  GROUP_SNV,
  GROUPS,
  DNA_OR_AA,
  COORDINATE_MODES,
  TABS,
} from '../constants/defs.json';

export const initialDataValues = {
  aggSequencesLocationGroupDate: [],
  aggSequencesGroupDate: [],
  dataAggLocationGroupDate: [],
  dataAggGroupDate: [],
  dataAggGroup: [],

  // Metadata filtering
  metadataCounts: {},

  dataAggLocationSnvDate: [],
  dataAggSnvDate: [],
  snvCooccurrence: [],

  countsPerLocation: {},
  countsPerLocationDate: {},
  validGroups: {},
  groupCounts: [],
};

let UIStoreInstance;
let configStoreInstance;
let snpDataStoreInstance;

export class DataStore {
  dataDate;
  numSequences;
  @observable numSequencesAfterAllFiltering;
  aggSequencesLocationGroupDate =
    initialDataValues.aggSequencesLocationGroupDate;
  aggSequencesGroupDate = initialDataValues.aggSequencesGroupDate;
  dataAggLocationGroupDate = initialDataValues.dataAggLocationGroupDate;
  dataAggGroupDate = initialDataValues.dataAggGroupDate;
  dataAggGroup = initialDataValues.dataAggGroup;
  @observable metadataCounts = initialDataValues.metadataCounts;

  dataAggLocationSnvDate = initialDataValues.dataAggLocationSnvDate;
  dataAggSnvDate = initialDataValues.dataAggSnvDate;
  snvCooccurrence = initialDataValues.snvCooccurrence;

  countsPerLocation = initialDataValues.countsPerLocation;
  countsPerLocationDate = initialDataValues.countsPerLocationDate;
  validGroups = initialDataValues.validGroups;
  groupCounts = initialDataValues.groupCounts;

  constructor() {}

  init() {
    this.dataDate = asyncDataStoreInstance.data.data_date;
    this.numSequences = asyncDataStoreInstance.data.num_sequences;

    UIStoreInstance = rootStoreInstance.UIStore;
    configStoreInstance = rootStoreInstance.configStore;
    snpDataStoreInstance = rootStoreInstance.snpDataStore;

    if (
      UIStoreInstance.activeTab === TABS.TAB_GROUP ||
      UIStoreInstance.activeTab === TABS.TAB_LOCATION
    ) {
      this.fetchData();
    }
  }

  @action
  fetchData() {
    UIStoreInstance.onCaseDataStateStarted();

    fetch(hostname + '/data', {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        group_key: toJS(configStoreInstance.groupKey),
        dna_or_aa: toJS(configStoreInstance.dnaOrAa),
        coordinate_mode: toJS(configStoreInstance.coordinateMode),
        coordinate_ranges: configStoreInstance.getCoordinateRanges(),
        selected_gene: toJS(configStoreInstance.selectedGene).name,
        selected_protein: toJS(configStoreInstance.selectedProtein).name,
        location_ids: getLocationIdsByNode(
          toJS(configStoreInstance.selectedLocationNodes)
        ),
        selected_metadata_fields: configStoreInstance.getSelectedMetadataFields(),
        ageRange: toJS(configStoreInstance.ageRange),
        low_count_filter: toJS(configStoreInstance.lowFreqFilterType),
        max_group_counts: parseInt(toJS(configStoreInstance.maxGroupCounts)),
        min_local_counts: parseInt(toJS(configStoreInstance.minLocalCounts)),
        min_global_counts: parseInt(toJS(configStoreInstance.minGlobalCounts)),
        start_date: toJS(configStoreInstance.startDate),
        end_date: toJS(configStoreInstance.endDate),
      }),
    })
      .then((res) => {
        if (!res.ok) {
          throw res;
        }
        return res.json();
      })
      .then((pkg) => {
        this.aggSequencesLocationGroupDate = pkg.aggSequencesLocationGroupDate;
        this.aggSequencesGroupDate = pkg.aggSequencesGroupDate;
        this.numSequencesAfterAllFiltering = pkg.numSequences;
        this.dataAggLocationGroupDate = pkg.dataAggLocationGroupDate;
        this.dataAggGroupDate = pkg.dataAggGroupDate;
        this.metadataCounts = pkg.metadataCounts;
        this.countsPerLocation = pkg.countsPerLocation;
        this.countsPerLocationDate = pkg.countsPerLocationDate;
        this.validGroups = pkg.validGroups;
        this.dataAggGroup = pkg.dataAggGroup;
        this.groupCounts = pkg.groupCounts;

        UIStoreInstance.onCaseDataStateFinished();

        if (configStoreInstance.groupKey === GROUP_SNV) {
          this.processSelectedSnvs();
        }
      })
      .catch((err) => {
        err
          .text()
          .then((errMsg) => {
            console.error(errMsg);
          })
          .finally(() => {
            UIStoreInstance.onCaseDataStateErr();
            console.error('Error fetching data');
          });
      });
  }

  @action
  processSelectedSnvs() {
    UIStoreInstance.onSnvDataStarted();
    const { snvColorMap } = snpDataStoreInstance;

    this.processCooccurrenceData();
    processSelectedSnvs(
      {
        selectedGroupIds: configStoreInstance.getSelectedGroupIds(),
        intToSnvMap: configStoreInstance.getIntToSnvMap(),
        dnaOrAa: toJS(configStoreInstance.dnaOrAa),
        countsPerLocation: this.countsPerLocation,
        validGroups: this.validGroups,
        aggSequencesLocationGroupDate: this.aggSequencesLocationGroupDate,
        aggSequencesGroupDate: this.aggSequencesGroupDate,
        // SNV data
        snvColorMap,
      },
      ({ dataAggLocationSnvDate, dataAggSnvDate }) => {
        this.dataAggLocationSnvDate = dataAggLocationSnvDate;
        this.dataAggSnvDate = dataAggSnvDate;
        UIStoreInstance.onSnvDataFinished();
      }
    );
  }

  @action
  processCooccurrenceData() {
    UIStoreInstance.onCooccurrenceDataStarted();

    const { snvColorMap } = snpDataStoreInstance;

    processCooccurrenceData(
      {
        selectedGroupIds: configStoreInstance.getSelectedGroupIds(),
        intToSnvMap: configStoreInstance.getIntToSnvMap(),
        dnaOrAa: toJS(configStoreInstance.dnaOrAa),
        aggSequencesGroupDate: this.aggSequencesGroupDate,
        // SNV data
        snvColorMap,
      },
      ({ snvCooccurrence }) => {
        this.snvCooccurrence = snvCooccurrence;
        UIStoreInstance.onCooccurrenceDataFinished();
      }
    );
  }

  @action
  downloadSelectedSequenceMetadata({ selectedFields, snvFormat }) {
    UIStoreInstance.onDownloadStarted();

    fetch(hostname + '/download_metadata', {
      method: 'POST',
      headers: {
        Accept: 'text/csv',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        location_ids: getLocationIdsByNode(
          toJS(configStoreInstance.selectedLocationNodes)
        ),
        selected_metadata_fields: configStoreInstance.getSelectedMetadataFields(),
        ageRange: toJS(configStoreInstance.ageRange),
        start_date: toJS(configStoreInstance.startDate),
        end_date: toJS(configStoreInstance.endDate),
        // Pass an array of only the fields that were selected
        selected_fields: Object.keys(selectedFields).filter(
          (field) => selectedFields[field]
        ),
        snv_format: snvFormat,
      }),
    })
      .then((res) => {
        if (!res.ok) {
          throw res;
        }
        return res.blob();
      })
      .then((blob) => {
        const url = URL.createObjectURL(blob);
        downloadBlobURL(url, 'selected_sequence_metadata.csv');
        UIStoreInstance.onDownloadFinished();
      })
      .catch((err) => {
        err
          .text()
          .then((errMsg) => {
            console.error(errMsg);
          })
          .finally(() => {
            UIStoreInstance.onDownloadErr();
            console.error('Download metadata failed');
          });
      });
  }

  downloadSelectedSNVs() {
    fetch(hostname + '/download_snvs', {
      method: 'POST',
      headers: {
        Accept: 'text/csv',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        group_key: toJS(configStoreInstance.groupKey),
        dna_or_aa: toJS(configStoreInstance.dnaOrAa),
        coordinate_mode: toJS(configStoreInstance.coordinateMode),
        coordinate_ranges: configStoreInstance.getCoordinateRanges(),
        selected_gene: toJS(configStoreInstance.selectedGene).name,
        selected_protein: toJS(configStoreInstance.selectedProtein).name,
        location_ids: getLocationIdsByNode(
          toJS(configStoreInstance.selectedLocationNodes)
        ),
        selected_metadata_fields: configStoreInstance.getSelectedMetadataFields(),
        ageRange: toJS(configStoreInstance.ageRange),
        start_date: toJS(configStoreInstance.startDate),
        end_date: toJS(configStoreInstance.endDate),
      }),
    })
      .then((res) => {
        if (!res.ok) {
          throw res;
        }
        return res.blob();
      })
      .then((blob) => {
        const url = URL.createObjectURL(blob);
        downloadBlobURL(url, 'selected_snvs.csv');
      })
      .catch((err) => {
        err
          .text()
          .then((errMsg) => {
            console.error(errMsg);
          })
          .finally(() => {
            UIStoreInstance.onDownloadErr();
            console.error('Download selected SNVs failed');
          });
      });
  }

  // TODO:
  // We should probably change this request to use form data
  // or query params, so that we can open the request in a new
  // window. Right now the gzipped file has to be loaded as a
  // blob on the front-end, and this will get unsustainable
  // if the user tries to download 100,000+ genomes
  @action
  downloadGenomes() {
    UIStoreInstance.onDownloadStarted();

    fetch(hostname + '/download_genomes', {
      method: 'POST',
      headers: {
        Accept: 'application/gzip',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        location_ids: getLocationIdsByNode(
          toJS(configStoreInstance.selectedLocationNodes)
        ),
        selected_metadata_fields: configStoreInstance.getSelectedMetadataFields(),
        ageRange: toJS(configStoreInstance.ageRange),
        start_date: toJS(configStoreInstance.startDate),
        end_date: toJS(configStoreInstance.endDate),
      }),
    })
      .then((res) => {
        if (!res.ok) {
          throw res;
        }
        return res.blob();
      })
      .then((blob) => {
        const url = URL.createObjectURL(blob);
        downloadBlobURL(url, 'genomes.fa.gz');
        UIStoreInstance.onDownloadFinished();
      })
      .catch((err) => {
        err
          .text()
          .then((errMsg) => {
            console.error(errMsg);
          })
          .finally(() => {
            UIStoreInstance.onDownloadErr();
            console.error('Error downloading genomes');
          });
      });
  }

  downloadAggSequences() {
    let csvString = `location,collection_date,${configStoreInstance.getGroupLabel()},count\n`;

    let intToSnvMap = configStoreInstance.getIntToSnvMap();
    this.aggSequencesLocationGroupDate.forEach((row) => {
      // Location data, date
      csvString += `"${row.location}","${intToISO(row.collection_date)}",`;

      // Group
      if (configStoreInstance.groupKey === GROUP_SNV) {
        csvString += `"${row.group_id
          .map((snvId) => intToSnvMap[snvId].snp_str)
          .join(';')}",`;
      } else {
        csvString += `"${row.group_id}",`;
      }

      // Counts
      csvString += `${row.counts}\n`;
    });

    const blob = new Blob([csvString]);
    const url = URL.createObjectURL(blob);

    downloadBlobURL(url, 'aggregate_sequences.csv');
  }

  downloadDataAggGroupDate() {
    // Write to a CSV string
    let csvString = `collection_date,${configStoreInstance.getGroupLabel()},count\n`;
    this.dataAggGroupDate.forEach((row) => {
      csvString += `${intToISO(parseInt(row.date))},${row.group},${
        row.counts
      }\n`;
    });

    const blob = new Blob([csvString]);
    const url = URL.createObjectURL(blob);

    downloadBlobURL(url, 'data_agg_group_date.csv');
  }

  downloadDataAggLocationGroupDate() {
    let locationData = JSON.parse(
      JSON.stringify(this.dataAggLocationGroupDate)
    );

    let csvString = `location,collection_date,${configStoreInstance.getGroupLabel()},${configStoreInstance.getGroupLabel()} Name,count\n`;

    locationData.forEach((row) => {
      csvString += `${row.location},${intToISO(parseInt(row.date))},${
        row.group
      },${row.group_name},${row.counts}\n`;
    });

    const blob = new Blob([csvString]);
    const url = URL.createObjectURL(blob);

    downloadBlobURL(url, 'data_agg_location_group_date.csv');
  }

  downloadSnvFrequencies() {
    let csvString = 'snv,';
    let fields = [];
    if (configStoreInstance.dnaOrAa === DNA_OR_AA.AA) {
      if (configStoreInstance.coordinateMode === COORDINATE_MODES.COORD_GENE) {
        csvString += 'gene,';
        fields.push('gene');
      } else if (
        configStoreInstance.coordinateMode === COORDINATE_MODES.COORD_PROTEIN
      ) {
        csvString += 'protein,';
        fields.push('protein');
      }
    }
    csvString += 'pos,ref,alt,counts\n';
    fields.push('pos', 'ref', 'alt');

    // Have to convert SNV string into integer, then into SNV object
    const snvToIntMap = configStoreInstance.getSnvToIntMap();
    const intToSnvMap = configStoreInstance.getIntToSnvMap();
    let snv;

    // groupCounts is a list of lists
    // [group: str, count: int, color: str, group_name: str]
    this.groupCounts
      .sort((a, b) => b[1] - a[1]) // Sort by counts, descending order
      .forEach((item) => {
        if (item[0] === GROUPS.OTHER_GROUP) {
          return;
        }
        // Add SNV name
        csvString += formatSnv(item[0], configStoreInstance.dnaOrAa) + ',';
        snv = intToSnvMap[snvToIntMap[item[0]]];
        // Add SNV fields
        if (item[0] === GROUPS.REFERENCE_GROUP) {
          csvString += fields.fill('').join(',');
        } else {
          csvString += fields.map((field) => snv[field]).join(',');
        }
        // Add counts
        csvString += `,${item[1]}\n`;
      });

    const blob = new Blob([csvString]);
    const url = URL.createObjectURL(blob);

    downloadBlobURL(url, 'snv_frequencies.csv');
  }

  downloadSnvCooccurrence() {
    let csvString = 'selected_snvs,cooccurs_with,';
    let fields = [];
    if (configStoreInstance.dnaOrAa === DNA_OR_AA.AA) {
      if (configStoreInstance.coordinateMode === COORDINATE_MODES.COORD_GENE) {
        csvString += 'gene,';
        fields.push('gene');
      } else if (
        configStoreInstance.coordinateMode === COORDINATE_MODES.COORD_PROTEIN
      ) {
        csvString += 'protein,';
        fields.push('protein');
      }
    }
    csvString += 'pos,ref,alt,count,selected_count,fraction\n';
    fields.push('pos', 'ref', 'alt');

    // Have to convert SNV string into integer, then into SNV object
    const snvToIntMap = configStoreInstance.getSnvToIntMap();
    const intToSnvMap = configStoreInstance.getIntToSnvMap();
    let snv;

    this.snvCooccurrence
      .sort((a, b) => {
        // Group by the SNV combination, then sort by counts, descending order
        if (a.combi < b.combi) {
          return -1;
        } else if (a.combi > b.combi) {
          return 1;
        } else {
          return b.count - a.count;
        }
      })
      .forEach((row) => {
        // Add combination and co-occuring SNV name
        csvString += `${row.combiName},${row.snvName},`;
        // Add SNV fields
        snv = intToSnvMap[snvToIntMap[row.snv]];
        csvString += fields.map((field) => snv[field]).join(',');

        // Add counts and fraction
        csvString += `,${row.count},${row.combiCount},${row.fraction}\n`;
      });

    const blob = new Blob([csvString]);
    const url = URL.createObjectURL(blob);

    downloadBlobURL(url, 'snv_cooccurrence.csv');
  }

  downloadCountryScoreData() {
    let jsonString = JSON.stringify(
      rootStoreInstance.globalSequencingDataStore.countryScoreData
    );
    const blob = new Blob([jsonString]);
    const url = URL.createObjectURL(blob);

    downloadBlobURL(url, 'global_sequencing_coverage.json');
  }
}

export default DataStore;
