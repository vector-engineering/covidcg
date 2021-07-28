import { action, toJS } from 'mobx';
import { hostname } from '../config';
import {
  processSelectedSnvs,
  processCooccurrenceData,
} from '../utils/snpDataWorkerWrapper';
import { downloadBlobURL } from '../utils/download';
import { intToISO } from '../utils/date';
import { getLocationIdsByNode } from '../utils/location';
import { asyncDataStoreInstance } from '../components/App';
import { rootStoreInstance } from './rootStore';
import {
  GROUP_SNV,
  GROUPS,
  DNA_OR_AA,
  COORDINATE_MODES,
  TABS,
} from '../constants/defs.json';

import {
  removeSubsetLocations,
  aggregateGroupDate,
  countGroups,
  getLocationCounts,
  expandSingleSnvData,
} from '../utils/data';

export class DataStore {
  dataDate;
  numSequences;
  numSequencesAfterAllFiltering;
  aggLocationGroupDate = [];
  aggGroupDate = [];
  aggLocationSingleSnvDate = [];

  aggLocationSelectedSnvsDate = [];
  aggSelectedSnvsDate = [];
  snvCooccurrence = [];

  countsPerLocationDateMap = new Map();
  cumulativeCountsPerLocationDateMap = new Map();
  countsPerLocationMap = {};
  groupCounts = [];

  constructor() {}

  init() {
    this.dataDate = asyncDataStoreInstance.data.data_date;
    this.numSequences = asyncDataStoreInstance.data.num_sequences;

    if (
      rootStoreInstance.UIStore.activeTab === TABS.TAB_COMPARE_GROUPS ||
      rootStoreInstance.UIStore.activeTab === TABS.TAB_COMPARE_LOCATIONS
    ) {
      this.fetchData();
    }
  }

  @action
  fetchData = () => {
    rootStoreInstance.UIStore.onCaseDataStateStarted();

    fetch(hostname + '/data', {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        group_key: toJS(rootStoreInstance.configStore.groupKey),
        dna_or_aa: toJS(rootStoreInstance.configStore.dnaOrAa),
        coordinate_mode: toJS(rootStoreInstance.configStore.coordinateMode),
        coordinate_ranges: rootStoreInstance.configStore.getCoordinateRanges(),
        selected_gene: toJS(rootStoreInstance.configStore.selectedGene).name,
        selected_protein: toJS(rootStoreInstance.configStore.selectedProtein)
          .name,
        location_ids: getLocationIdsByNode(
          toJS(rootStoreInstance.configStore.selectedLocationNodes)
        ),
        selected_metadata_fields:
          rootStoreInstance.configStore.getSelectedMetadataFields(),
        ageRange: toJS(rootStoreInstance.configStore.ageRange),
        start_date: toJS(rootStoreInstance.configStore.startDate),
        end_date: toJS(rootStoreInstance.configStore.endDate),
      }),
    })
      .then((res) => {
        if (!res.ok) {
          throw res;
        }
        return res.json();
      })
      .then((res) => {
        this.aggLocationGroupDate = res;

        // console.log(this.aggLocationGroupDate);

        // Create copy of the data with subset locations removed
        this.aggSequencesUniqueLocationGroupDate = removeSubsetLocations({
          aggLocationGroupDate: this.aggLocationGroupDate,
          selectedLocationNodes: toJS(
            rootStoreInstance.configStore.selectedLocationNodes
          ),
        });

        // Count all sequences
        this.numSequencesAfterAllFiltering =
          this.aggSequencesUniqueLocationGroupDate.reduce(
            (accumulator, record) => accumulator + record.counts,
            0
          );

        // Count groups
        this.groupCounts = countGroups({
          aggLocationGroupDate: this.aggSequencesUniqueLocationGroupDate,
          groupKey: toJS(rootStoreInstance.configStore.groupKey),
        });

        // Collapse low frequency groups
        // Identify groups to collapse into the "Other" group
        // this.validGroups = getValidGroups({
        //   aggSequencesGroup: this.aggSequencesGroup,
        //   lowFreqFilterType: toJS(rootStoreInstance.configStore.lowFreqFilterType),
        //   lowFreqFilterParams: {
        //     maxGroupCounts: toJS(rootStoreInstance.configStore.maxGroupCounts),
        //     minLocalCounts: toJS(rootStoreInstance.configStore.minLocalCounts),
        //   },
        // });
        // Make a copy of the data which has collapsed data

        // console.log(this.groupCounts.sort((a, b) => a.group_id - b.group_id));

        ({ aggGroupDate: this.aggGroupDate, aggGroup: this.aggSequencesGroup } =
          aggregateGroupDate({
            aggSequencesUniqueLocationGroupDate:
              this.aggSequencesUniqueLocationGroupDate,
            groupKey: rootStoreInstance.configStore.groupKey,
          }));
        // console.log(this.aggGroupDate);
        // console.log(this.aggSequencesGroup);

        ({
          countsPerLocationDateMap: this.countsPerLocationDateMap,
          cumulativeCountsPerLocationDateMap:
            this.cumulativeCountsPerLocationDateMap,
          countsPerLocationMap: this.countsPerLocationMap,
        } = getLocationCounts({
          aggLocationGroupDate: this.aggLocationGroupDate,
        }));

        // Expand aggLocationGroupDate into
        // single SNV data
        // i.e., transform data where each row represents
        // a co-occurring SNV, into data where each row represents
        // individual SNVs
        if (rootStoreInstance.configStore.groupKey === GROUP_SNV) {
          this.aggLocationSingleSnvDate = expandSingleSnvData({
            aggLocationGroupDate: this.aggLocationGroupDate,
          });
          // console.log(this.aggLocationSingleSnvDate);
        }

        rootStoreInstance.UIStore.onCaseDataStateFinished();

        if (rootStoreInstance.configStore.groupKey === GROUP_SNV) {
          this.processSelectedSnvs();
        }
      })
      .catch((err) => {
        let prefix = 'Error fetching data';
        if (!(typeof err.text === 'function')) {
          console.error(prefix, err);
        } else {
          err.text().then((errMsg) => {
            console.error(prefix, errMsg);
          });
        }
        rootStoreInstance.UIStore.onCaseDataStateErr();
      });
  };

  @action
  processSelectedSnvs = () => {
    rootStoreInstance.UIStore.onSnvDataStarted();
    const { snvColorMap } = rootStoreInstance.snpDataStore;

    this.processCooccurrenceData();
    processSelectedSnvs(
      {
        selectedGroupIds: rootStoreInstance.configStore.getSelectedGroupIds(),
        intToSnvMap: rootStoreInstance.configStore.getIntToSnvMap(),
        dnaOrAa: toJS(rootStoreInstance.configStore.dnaOrAa),
        countsPerLocationMap: this.countsPerLocationMap,
        // validGroups: this.validGroups,
        aggLocationGroupDate: this.aggLocationGroupDate,
        aggGroupDate: this.aggGroupDate,
        // SNV data
        snvColorMap,
      },
      ({ aggLocationSelectedSnvsDate, aggSelectedSnvsDate }) => {
        this.aggLocationSelectedSnvsDate = aggLocationSelectedSnvsDate;
        this.aggSelectedSnvsDate = aggSelectedSnvsDate;
        rootStoreInstance.UIStore.onSnvDataFinished();
      }
    );
  };

  @action
  processCooccurrenceData = () => {
    rootStoreInstance.UIStore.onCooccurrenceDataStarted();

    const { snvColorMap } = rootStoreInstance.snpDataStore;

    processCooccurrenceData(
      {
        selectedGroupIds: rootStoreInstance.configStore.getSelectedGroupIds(),
        intToSnvMap: rootStoreInstance.configStore.getIntToSnvMap(),
        dnaOrAa: toJS(rootStoreInstance.configStore.dnaOrAa),
        aggSequencesGroup: this.aggSequencesGroup,
        // SNV data
        snvColorMap,
      },
      ({ snvCooccurrence }) => {
        this.snvCooccurrence = snvCooccurrence;
        rootStoreInstance.UIStore.onCooccurrenceDataFinished();
      }
    );
  };

  @action
  downloadSelectedSequenceMetadata = ({ selectedFields, snvFormat }) => {
    rootStoreInstance.UIStore.onDownloadStarted();

    fetch(hostname + '/download_metadata', {
      method: 'POST',
      headers: {
        Accept: 'text/csv',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        location_ids: getLocationIdsByNode(
          toJS(rootStoreInstance.configStore.selectedLocationNodes)
        ),
        selected_metadata_fields:
          rootStoreInstance.configStore.getSelectedMetadataFields(),
        ageRange: toJS(rootStoreInstance.configStore.ageRange),
        start_date: toJS(rootStoreInstance.configStore.startDate),
        end_date: toJS(rootStoreInstance.configStore.endDate),
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
        rootStoreInstance.UIStore.onDownloadFinished();
      })
      .catch((err) => {
        let prefix = 'Error downloading metadata';
        if (!(typeof err.text === 'function')) {
          console.error(prefix, err);
        } else {
          err.text().then((errMsg) => {
            console.error(prefix, errMsg);
          });
        }
        rootStoreInstance.UIStore.onDownloadErr();
      });
  };

  @action
  downloadSelectedSNVs = () => {
    fetch(hostname + '/download_snvs', {
      method: 'POST',
      headers: {
        Accept: 'text/csv',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        group_key: toJS(rootStoreInstance.configStore.groupKey),
        dna_or_aa: toJS(rootStoreInstance.configStore.dnaOrAa),
        coordinate_mode: toJS(rootStoreInstance.configStore.coordinateMode),
        coordinate_ranges: rootStoreInstance.configStore.getCoordinateRanges(),
        selected_gene: toJS(rootStoreInstance.configStore.selectedGene).name,
        selected_protein: toJS(rootStoreInstance.configStore.selectedProtein)
          .name,
        location_ids: getLocationIdsByNode(
          toJS(rootStoreInstance.configStore.selectedLocationNodes)
        ),
        selected_metadata_fields:
          rootStoreInstance.configStore.getSelectedMetadataFields(),
        ageRange: toJS(rootStoreInstance.configStore.ageRange),
        start_date: toJS(rootStoreInstance.configStore.startDate),
        end_date: toJS(rootStoreInstance.configStore.endDate),
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
        let prefix = 'Download selected SNVs failed';
        if (!(typeof err.text === 'function')) {
          console.error(prefix, err);
        } else {
          err.text().then((errMsg) => {
            console.error(prefix, errMsg);
          });
        }
        rootStoreInstance.UIStore.onDownloadErr();
      });
  };

  // TODO:
  // We should probably change this request to use form data
  // or query params, so that we can open the request in a new
  // window. Right now the gzipped file has to be loaded as a
  // blob on the front-end, and this will get unsustainable
  // if the user tries to download 100,000+ genomes
  @action
  downloadGenomes = () => {
    rootStoreInstance.UIStore.onDownloadStarted();

    fetch(hostname + '/download_genomes', {
      method: 'POST',
      headers: {
        Accept: 'application/gzip',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        location_ids: getLocationIdsByNode(
          toJS(rootStoreInstance.configStore.selectedLocationNodes)
        ),
        selected_metadata_fields:
          rootStoreInstance.configStore.getSelectedMetadataFields(),
        ageRange: toJS(rootStoreInstance.configStore.ageRange),
        start_date: toJS(rootStoreInstance.configStore.startDate),
        end_date: toJS(rootStoreInstance.configStore.endDate),
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
        rootStoreInstance.UIStore.onDownloadFinished();
      })
      .catch((err) => {
        let prefix = 'Error downloading genomes';
        if (!(typeof err.text === 'function')) {
          console.error(prefix, err);
        } else {
          err.text().then((errMsg) => {
            console.error(prefix, errMsg);
          });
        }
        rootStoreInstance.UIStore.onDownloadErr();
      });
  };

  downloadAggSequences() {
    let csvString = `location,collection_date,${rootStoreInstance.configStore.getGroupLabel()},count\n`;

    let intToSnvMap = rootStoreInstance.configStore.getIntToSnvMap();
    this.aggLocationGroupDate.forEach((row) => {
      // Location data, date
      csvString += `"${row.location}","${intToISO(row.collection_date)}",`;

      // Group
      if (rootStoreInstance.configStore.groupKey === GROUP_SNV) {
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

  downloadAggGroupDate() {
    // Write to a CSV string
    let csvString = `collection_date,"${rootStoreInstance.configStore.getGroupLabel()}",count\n`;
    this.aggGroupDate.forEach((row) => {
      let groupName;
      // If we're in SNV mode, then we have to map SNV IDs back
      // to a co-occurring SNV string representation
      if (rootStoreInstance.configStore.groupKey === GROUP_SNV) {
        groupName = row.group_id
          .map(
            (snvId) =>
              rootStoreInstance.snpDataStore.intToSnv(
                rootStoreInstance.configStore.dnaOrAa,
                rootStoreInstance.configStore.coordinateMode,
                snvId
              ).name
          )
          .join(';');
      } else {
        groupName = row.group_id;
      }

      csvString += `${intToISO(row.collection_date)},"${groupName}",${
        row.counts
      }\n`;
    });

    const blob = new Blob([csvString]);
    const url = URL.createObjectURL(blob);

    downloadBlobURL(url, 'data_agg_group_date.csv');
  }

  downloadAggLocationGroupDate() {
    let locationData;
    if (rootStoreInstance.configStore.groupKey === GROUP_SNV) {
      locationData = toJS(this.aggLocationSingleSnvDate);
      // Get SNV data
      locationData.forEach((record) => {
        let snv = rootStoreInstance.snpDataStore.intToSnv(
          rootStoreInstance.configStore.dnaOrAa,
          rootStoreInstance.configStore.coordinateMode,
          record.group_id
        );
        record.group_name = snv.name;
        record.group = snv.snp_str;
      });
    } else {
      locationData = toJS(this.aggLocationGroupDate);
      locationData.forEach((record) => {
        record.group_name = record.group_id;
        record.group = record.group_id;
      });
    }
    let csvString = `location,collection_date,${rootStoreInstance.configStore.getGroupLabel()},${rootStoreInstance.configStore.getGroupLabel()} Name,count\n`;
    locationData.forEach((row) => {
      csvString += `${row.location},${intToISO(row.collection_date)},${
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
    if (rootStoreInstance.configStore.dnaOrAa === DNA_OR_AA.AA) {
      if (
        rootStoreInstance.configStore.coordinateMode ===
        COORDINATE_MODES.COORD_GENE
      ) {
        csvString += 'gene,';
        fields.push('gene');
      } else if (
        rootStoreInstance.configStore.coordinateMode ===
        COORDINATE_MODES.COORD_PROTEIN
      ) {
        csvString += 'protein,';
        fields.push('protein');
      }
    }
    csvString += 'pos,ref,alt,counts\n';
    fields.push('pos', 'ref', 'alt');

    let snv;

    // groupCounts is a list of records
    // [{group_id: int, group: str, counts: int, color: str, group_name: str}]
    this.groupCounts
      .sort((a, b) => b.counts - a.counts) // Sort by counts, descending order
      .forEach((record) => {
        snv = rootStoreInstance.snpDataStore.intToSnv(
          rootStoreInstance.configStore.dnaOrAa,
          rootStoreInstance.configStore.coordinateMode,
          record.group_id
        );

        // Add SNV fields
        if (snv.snp_str === GROUPS.REFERENCE_GROUP) {
          csvString += snv.snp_str + ',';
          csvString += fields.slice().fill('').join(',');
        } else {
          csvString +=
            snv.name + ',' + fields.map((field) => snv[field]).join(',');
        }
        // Add counts
        csvString += `,${record.counts}\n`;
      });
    // console.log(csvString);

    const blob = new Blob([csvString]);
    const url = URL.createObjectURL(blob);

    downloadBlobURL(url, 'snv_frequencies.csv');
  }

  downloadSnvCooccurrence() {
    let csvString = 'selected_snvs,cooccurs_with,';
    let fields = [];
    if (rootStoreInstance.configStore.dnaOrAa === DNA_OR_AA.AA) {
      if (
        rootStoreInstance.configStore.coordinateMode ===
        COORDINATE_MODES.COORD_GENE
      ) {
        csvString += 'gene,';
        fields.push('gene');
      } else if (
        rootStoreInstance.configStore.coordinateMode ===
        COORDINATE_MODES.COORD_PROTEIN
      ) {
        csvString += 'protein,';
        fields.push('protein');
      }
    }
    csvString += 'pos,ref,alt,count,selected_count,fraction\n';
    fields.push('pos', 'ref', 'alt');

    // Have to convert SNV string into integer, then into SNV object
    const snvToIntMap = rootStoreInstance.configStore.getSnvToIntMap();
    const intToSnvMap = rootStoreInstance.configStore.getIntToSnvMap();
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
