import { observable, action, toJS } from 'mobx';
import { hostname } from '../config';
import {
  processSelectedSnvs,
  processCooccurrenceData,
} from '../utils/snpDataWorkerWrapper';
import { downloadBlobURL } from '../utils/download';
import { aggregate } from '../utils/transform';
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
} from '../constants/defs.json';

export const initialDataValues = {
  aggSequences: [],
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

export class DataStore {
  // References to store instances
  UIStoreInstance;
  configStoreInstance;
  snpDataStoreInstance;

  dataDate;
  numSequences;
  @observable numSequencesAfterAllFiltering;
  aggSequences = initialDataValues.aggSequences;
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

    this.UIStoreInstance = rootStoreInstance.UIStore;
    this.configStoreInstance = rootStoreInstance.configStore;
    this.snpDataStoreInstance = rootStoreInstance.snpDataStore;
    this.metadataStoreInstance = rootStoreInstance.metadataStore;
    this.locationStoreInstance = rootStoreInstance.locationDataStore;

    this.fetchData();
  }

  @action
  emptyCaseData() {
    this.UIStoreInstance.onCaseDataStateStarted();

    this.aggSequences = initialDataValues.aggSequences;
    this.dataAggLocationGroupDate = initialDataValues.dataAggLocationGroupDate;
    this.dataAggGroupDate = initialDataValues.dataAggGroupDate;
    this.dataAggGroup = initialDataValues.dataAggGroup;

    this.dataAggLocationSnvDate = initialDataValues.dataAggLocationSnvDate;
    this.dataAggSnvDate = initialDataValues.dataAggSnvDate;
    this.snvCooccurrence = initialDataValues.snvCooccurrence;

    this.countsPerLocation = initialDataValues.countsPerLocation;
    this.countsPerLocationDate = initialDataValues.countsPerLocationDate;

    this.UIStoreInstance.onCaseDataStateFinished();
  }

  @action
  async fetchData() {
    this.UIStoreInstance.onCaseDataStateStarted();

    const res = await fetch(hostname + '/data', {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        group_key: toJS(this.configStoreInstance.groupKey),
        dna_or_aa: toJS(this.configStoreInstance.dnaOrAa),
        coordinate_mode: toJS(this.configStoreInstance.coordinateMode),
        coordinate_ranges: this.configStoreInstance.getCoordinateRanges(),
        selected_gene: toJS(this.configStoreInstance.selectedGene).name,
        selected_protein: toJS(this.configStoreInstance.selectedProtein).name,
        location_ids: getLocationIdsByNode(
          toJS(this.configStoreInstance.selectedLocationNodes)
        ),
        selected_metadata_fields: this.configStoreInstance.getSelectedMetadataFields(),
        ageRange: toJS(this.configStoreInstance.ageRange),
        low_count_filter: toJS(this.configStoreInstance.lowFreqFilterType),
        max_group_counts: parseInt(
          toJS(this.configStoreInstance.maxGroupCounts)
        ),
        min_local_counts: parseInt(
          toJS(this.configStoreInstance.minLocalCounts)
        ),
        min_global_counts: parseInt(
          toJS(this.configStoreInstance.minGlobalCounts)
        ),
        start_date: toJS(this.configStoreInstance.startDate),
        end_date: toJS(this.configStoreInstance.endDate),
      }),
    });
    const pkg = await res.json();
    // console.log(pkg);

    this.aggSequences = pkg.aggSequences;
    this.numSequencesAfterAllFiltering = pkg.numSequences;
    this.dataAggLocationGroupDate = pkg.dataAggLocationGroupDate;
    this.dataAggGroupDate = pkg.dataAggGroupDate;
    this.metadataCounts = pkg.metadataCounts;
    this.countsPerLocation = pkg.countsPerLocation;
    this.countsPerLocationDate = pkg.countsPerLocationDate;
    this.validGroups = pkg.validGroups;
    this.dataAggGroup = pkg.dataAggGroup;
    this.groupCounts = pkg.groupCounts;

    this.UIStoreInstance.onCaseDataStateFinished();

    if (this.configStoreInstance.groupKey === GROUP_SNV) {
      this.processSelectedSnvs();
    }
  }

  @action
  processSelectedSnvs() {
    this.UIStoreInstance.onSnvDataStarted();
    const { snvColorMap } = this.snpDataStoreInstance;

    this.processCooccurrenceData();
    processSelectedSnvs(
      {
        selectedGroupIds: this.configStoreInstance.getSelectedGroupIds(),
        intToSnvMap: this.configStoreInstance.getIntToSnvMap(),
        dnaOrAa: toJS(this.configStoreInstance.dnaOrAa),
        countsPerLocation: this.countsPerLocation,
        validGroups: this.validGroups,
        aggSequences: this.aggSequences,
        // SNV data
        snvColorMap,
      },
      ({ dataAggLocationSnvDate, dataAggSnvDate }) => {
        this.dataAggLocationSnvDate = dataAggLocationSnvDate;
        this.dataAggSnvDate = dataAggSnvDate;
        this.UIStoreInstance.onSnvDataFinished();
      }
    );
  }

  @action
  processCooccurrenceData() {
    this.UIStoreInstance.onCooccurrenceDataStarted();

    const { snvColorMap } = this.snpDataStoreInstance;

    processCooccurrenceData(
      {
        selectedGroupIds: this.configStoreInstance.getSelectedGroupIds(),
        intToSnvMap: this.configStoreInstance.getIntToSnvMap(),
        dnaOrAa: toJS(this.configStoreInstance.dnaOrAa),
        aggSequences: this.aggSequences,
        // SNV data
        snvColorMap,
      },
      ({ snvCooccurrence }) => {
        this.snvCooccurrence = snvCooccurrence;
        this.UIStoreInstance.onCooccurrenceDataFinished();
      }
    );
  }

  @action
  async downloadSelectedSequenceMetadata({ selectedFields, snvFormat }) {
    this.UIStoreInstance.onDownloadStarted();

    const res = await fetch(hostname + '/download_metadata', {
      method: 'POST',
      headers: {
        Accept: 'text/csv',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        location_ids: getLocationIdsByNode(
          toJS(this.configStoreInstance.selectedLocationNodes)
        ),
        selected_metadata_fields: this.configStoreInstance.getSelectedMetadataFields(),
        ageRange: toJS(this.configStoreInstance.ageRange),
        start_date: toJS(this.configStoreInstance.startDate),
        end_date: toJS(this.configStoreInstance.endDate),
        // Pass an array of only the fields that were selected
        selected_fields: Object.keys(selectedFields).filter(
          (field) => selectedFields[field]
        ),
        snv_format: snvFormat,
      }),
    });
    const blob = await res.blob();
    const url = URL.createObjectURL(blob);
    downloadBlobURL(url, 'selected_sequence_metadata.csv');

    this.UIStoreInstance.onDownloadFinished();
  }

  async downloadSelectedSNVs() {
    const res = await fetch(hostname + '/download_snvs', {
      method: 'POST',
      headers: {
        Accept: 'text/csv',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        group_key: toJS(this.configStoreInstance.groupKey),
        dna_or_aa: toJS(this.configStoreInstance.dnaOrAa),
        coordinate_mode: toJS(this.configStoreInstance.coordinateMode),
        coordinate_ranges: this.configStoreInstance.getCoordinateRanges(),
        selected_gene: toJS(this.configStoreInstance.selectedGene).name,
        selected_protein: toJS(this.configStoreInstance.selectedProtein).name,
        location_ids: getLocationIdsByNode(
          toJS(this.configStoreInstance.selectedLocationNodes)
        ),
        selected_metadata_fields: this.configStoreInstance.getSelectedMetadataFields(),
        ageRange: toJS(this.configStoreInstance.ageRange),
        start_date: toJS(this.configStoreInstance.startDate),
        end_date: toJS(this.configStoreInstance.endDate),
      }),
    });
    const blob = await res.blob();
    const url = URL.createObjectURL(blob);
    downloadBlobURL(url, 'selected_snvs.csv');
  }

  // TODO:
  // We should probably change this request to use form data
  // or query params, so that we can open the request in a new
  // window. Right now the gzipped file has to be loaded as a
  // blob on the front-end, and this will get unsustainable
  // if the user tries to download 100,000+ genomes
  @action
  async downloadGenomes() {
    this.UIStoreInstance.onDownloadStarted();

    const res = await fetch(hostname + '/download_genomes', {
      method: 'POST',
      headers: {
        Accept: 'application/gzip',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        location_ids: getLocationIdsByNode(
          toJS(this.configStoreInstance.selectedLocationNodes)
        ),
        selected_metadata_fields: this.configStoreInstance.getSelectedMetadataFields(),
        ageRange: toJS(this.configStoreInstance.ageRange),
        start_date: toJS(this.configStoreInstance.startDate),
        end_date: toJS(this.configStoreInstance.endDate),
      }),
    });
    const blob = await res.blob();
    const url = URL.createObjectURL(blob);
    downloadBlobURL(url, 'genomes.fa.gz');

    this.UIStoreInstance.onDownloadFinished();
  }

  downloadAggSequences() {
    let csvString = `location,collection_date,${this.configStoreInstance.getGroupLabel()},count\n`;

    let intToSnvMap = this.configStoreInstance.getIntToSnvMap();
    this.aggSequences.forEach((row) => {
      // Location data, date
      csvString += `"${row.location}","${intToISO(row.collection_date)}",`;

      // Group
      if (this.configStoreInstance.groupKey === GROUP_SNV) {
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
    let csvString = `collection_date,${this.configStoreInstance.getGroupLabel()},count\n`;
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
    // console.log(this.dataAggLocationGroupDate);

    let locationData = JSON.parse(
      JSON.stringify(this.dataAggLocationGroupDate)
    );

    // Filter by date
    // if (
    //   this.configStoreInstance.dateRange[0] != -1 ||
    //   this.configStoreInstance.dateRange[1] != -1
    // ) {
    //   locationData = locationData.filter((row) => {
    //     return (
    //       (this.configStoreInstance.dateRange[0] == -1 ||
    //         row.date > this.configStoreInstance.dateRange[0]) &&
    //       (this.configStoreInstance.dateRange[1] == -1 ||
    //         row.date < this.configStoreInstance.dateRange[1])
    //     );
    //   });
    // }

    locationData = aggregate({
      data: locationData,
      groupby: ['location', 'date', 'group', 'group_name'],
      fields: ['counts', 'location_counts'],
      ops: ['sum', 'max'],
      as: ['counts', 'location_counts'],
    });

    // Manually join the countsPerLocationDate to locationData
    locationData.forEach((row) => {
      row.location_date_count = this.countsPerLocationDate[row.location][
        row.date.toString()
      ];
    });

    let csvString = `location,collection_date,${this.configStoreInstance.getGroupLabel()},count,location_date_count\n`;

    locationData.forEach((row) => {
      csvString += `${row.location},${intToISO(parseInt(row.date))},${
        row.group
      },${row.counts},${row.location_date_count}\n`;
    });

    const blob = new Blob([csvString]);
    const url = URL.createObjectURL(blob);

    downloadBlobURL(url, 'data_agg_location_group_date.csv');
  }

  downloadSnvFrequencies() {
    let csvString = 'snv,';
    let fields = [];
    if (this.configStoreInstance.dnaOrAa === DNA_OR_AA.AA) {
      if (
        this.configStoreInstance.coordinateMode === COORDINATE_MODES.COORD_GENE
      ) {
        csvString += 'gene,';
        fields.push('gene');
      } else if (
        this.configStoreInstance.coordinateMode ===
        COORDINATE_MODES.COORD_PROTEIN
      ) {
        csvString += 'protein,';
        fields.push('protein');
      }
    }
    csvString += 'pos,ref,alt,counts\n';
    fields.push('pos', 'ref', 'alt');

    // Have to convert SNV string into integer, then into SNV object
    const snvToIntMap = this.configStoreInstance.getSnvToIntMap();
    const intToSnvMap = this.configStoreInstance.getIntToSnvMap();
    let snv;

    this.groupCounts
      .sort((a, b) => b[1] - a[1]) // Sort by counts, descending order
      .forEach((item) => {
        if (item[0] === GROUPS.OTHER_GROUP) {
          return;
        }
        // Add SNV name
        csvString += formatSnv(item[0], this.configStoreInstance.dnaOrAa) + ',';
        snv = intToSnvMap[snvToIntMap[item[0]]];
        // Add SNV fields
        csvString += fields.map((field) => snv[field]).join(',');
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
    if (this.configStoreInstance.dnaOrAa === DNA_OR_AA.AA) {
      if (
        this.configStoreInstance.coordinateMode === COORDINATE_MODES.COORD_GENE
      ) {
        csvString += 'gene,';
        fields.push('gene');
      } else if (
        this.configStoreInstance.coordinateMode ===
        COORDINATE_MODES.COORD_PROTEIN
      ) {
        csvString += 'protein,';
        fields.push('protein');
      }
    }
    csvString += 'pos,ref,alt,count,selected_count,fraction\n';
    fields.push('pos', 'ref', 'alt');

    // Have to convert SNV string into integer, then into SNV object
    const snvToIntMap = this.configStoreInstance.getSnvToIntMap();
    const intToSnvMap = this.configStoreInstance.getIntToSnvMap();
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
    let jsonString = JSON.stringify(asyncDataStoreInstance.data.country_score);
    const blob = new Blob([jsonString]);
    const url = URL.createObjectURL(blob);

    downloadBlobURL(url, 'global_sequencing_coverage.json');
  }
}

export default DataStore;
