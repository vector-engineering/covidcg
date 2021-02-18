import { observable, action, toJS } from 'mobx';
import {
  processSelectedSnvs,
  processCooccurrenceData,
} from '../utils/snpDataWorkerWrapper';
import { downloadBlobURL, generateSelectionString } from '../utils/download';
import { aggregate } from '../utils/transform';
import { intToISO } from '../utils/date';
import { getLocationIdsByNode } from '../utils/location';
import { asyncDataStoreInstance } from '../components/App';
import { rootStoreInstance } from './rootStore';

export const initialDataValues = {
  filteredCaseData: [],
  dataAggLocationGroupDate: [],
  dataAggGroupDate: [],
  dataAggGroup: [],

  // Metadata filtering
  numSequencesBeforeMetadataFiltering: 0,
  metadataCounts: {},
  metadataCountsAfterFiltering: {},

  dataAggLocationSnvDate: [],
  dataAggSnvDate: [],
  snvCooccurrence: [],

  countsPerLocation: {},
  countsPerLocationDate: {},
  validGroups: {},
  countsPerGroup: {},
  countsPerGroupDateFiltered: [],
};

export class DataStore {
  // References to store instances
  UIStoreInstance;
  configStoreInstance;
  snpDataStoreInstance;

  dataDate;
  numSequences;
  @observable numSequencesAfterAllFiltering;
  filteredCaseData = initialDataValues.filteredCaseData;
  dataAggLocationGroupDate = initialDataValues.dataAggLocationGroupDate;
  dataAggGroupDate = initialDataValues.dataAggGroupDate;
  dataAggGroup = initialDataValues.dataAggGroup;
  @observable numSequencesBeforeMetadataFiltering =
    initialDataValues.numSequencesBeforeMetadataFiltering;
  @observable metadataCounts = initialDataValues.metadataCounts;
  @observable metadataCountsAfterFiltering =
    initialDataValues.metadataCountsAfterFiltering;

  dataAggLocationSnvDate = initialDataValues.dataAggLocationSnvDate;
  dataAggSnvDate = initialDataValues.dataAggSnvDate;
  snvCooccurrence = initialDataValues.snvCooccurrence;

  countsPerLocation = initialDataValues.countsPerLocation;
  countsPerLocationDate = initialDataValues.countsPerLocationDate;
  validGroups = initialDataValues.validGroups;
  countsPerGroup = initialDataValues.countsPerGroup;
  countsPerGroupDateFiltered = initialDataValues.countsPerGroupDateFiltered;

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

    const selectedMetadataFields = toJS(
      this.configStoreInstance.selectedMetadataFields
    );
    Object.keys(selectedMetadataFields).forEach((metadataField) => {
      selectedMetadataFields[metadataField] = selectedMetadataFields[
        metadataField
      ].map((item) => {
        return parseInt(item.value);
      });
    });

    const res = await fetch('http://localhost:5000/data', {
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
        selected_gene: toJS(this.configStoreInstance.selectedGene),
        selected_protein: toJS(this.configStoreInstance.selectedProtein),
        location_ids: getLocationIdsByNode(
          toJS(this.configStoreInstance.selectedLocationNodes)
        ),
        selected_metadata_fields: selectedMetadataFields,
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

    this.filteredCaseData = pkg.filteredCaseData;
    this.numSequencesAfterAllFiltering = pkg.filteredCaseData.length;
    this.dataAggLocationGroupDate = pkg.dataAggLocationGroupDate;
    this.dataAggGroupDate = pkg.dataAggGroupDate;
    this.metadataCounts = pkg.metadataCounts;
    this.metadataCountsAfterFiltering = pkg.metadataCountsAfterFiltering;
    this.numSequencesBeforeMetadataFiltering =
      pkg.numSequencesBeforeMetadataFiltering;
    this.countsPerLocation = pkg.countsPerLocation;
    this.countsPerLocationDate = pkg.countsPerLocationDate;
    this.validGroups = pkg.validGroups;
    this.countsPerGroup = pkg.countsPerGroup;
    this.dataAggGroup = pkg.dataAggGroup;
    this.countsPerGroupDateFiltered = pkg.countsPerGroupDateFiltered;

    this.UIStoreInstance.onCaseDataStateFinished();
  }

  @action
  processSelectedSnvs() {
    this.UIStoreInstance.onSnvDataStarted();
    const {
      intToDnaSnvMap,
      intToGeneAaSnvMap,
      intToProteinAaSnvMap,
      dnaSnvMap,
      geneAaSnvMap,
      proteinAaSnvMap,
      snvColorMap,
    } = this.snpDataStoreInstance;

    this.processCooccurrenceData();
    processSelectedSnvs(
      {
        dnaOrAa: toJS(this.configStoreInstance.dnaOrAa),
        coordinateMode: toJS(this.configStoreInstance.coordinateMode),
        selectedLocationNodes: toJS(
          this.configStoreInstance.selectedLocationNodes
        ),
        countsPerLocation: this.countsPerLocation,
        filteredCaseData: this.filteredCaseData,
        selectedGroups: toJS(this.configStoreInstance.selectedGroups).map(
          (item) => item.group
        ),
        validGroups: this.validGroups,

        // SNV data
        intToDnaSnvMap,
        intToGeneAaSnvMap,
        intToProteinAaSnvMap,
        dnaSnvMap,
        geneAaSnvMap,
        proteinAaSnvMap,
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

    const {
      intToDnaSnvMap,
      intToGeneAaSnvMap,
      intToProteinAaSnvMap,
      dnaSnvMap,
      geneAaSnvMap,
      proteinAaSnvMap,
      snvColorMap,
    } = this.snpDataStoreInstance;

    processCooccurrenceData(
      {
        dnaOrAa: toJS(this.configStoreInstance.dnaOrAa),
        coordinateMode: toJS(this.configStoreInstance.coordinateMode),
        filteredCaseData: this.filteredCaseData,
        selectedGroups: toJS(this.configStoreInstance.selectedGroups).map(
          (item) => item.group
        ),
        dateRange: toJS(this.configStoreInstance.dateRange),

        // SNV data
        intToDnaSnvMap,
        intToGeneAaSnvMap,
        intToProteinAaSnvMap,
        dnaSnvMap,
        geneAaSnvMap,
        proteinAaSnvMap,
        snvColorMap,
      },
      ({ snvCooccurrence }) => {
        this.snvCooccurrence = snvCooccurrence;
        this.UIStoreInstance.onCooccurrenceDataFinished();
      }
    );
  }

  @action
  downloadSelectedSequenceMetadata() {
    let csvString = [];
    // Take the first sequence to get the keys/fields for this file
    const keys = Object.keys(this.rawCaseData[0]);
    // Special case - we're going to expand the location id into
    // 4 fields (region, country, division, location)
    // so add these to the end
    csvString.push(
      keys.concat(['region', 'country', 'division', 'location']).join(',')
    );

    const metadataFields = Object.keys(this.metadataStoreInstance.metadataMap);

    let vals = [];
    let metadataVal;
    this.filteredCaseData.forEach((row) => {
      //  vals = keys.map(key => row[key]);
      vals = [];
      keys.forEach((key) => {
        // Special cases for keys
        if (key === 'collection_date') {
          vals.push('"' + intToISO(row[key]) + '"');
        } else if (key === 'age_start' || key === 'age_end') {
          vals.push(row[key] === null ? '' : row[key]);
        } else if (metadataFields.includes(key)) {
          metadataVal = this.metadataStoreInstance.getMetadataValueFromId(
            key,
            row[key]
          );
          vals.push('"' + (metadataVal === undefined ? '' : metadataVal) + '"');
        } else if (key === 'dna_snp_str') {
          vals.push(
            '"' +
              row[key]
                .map((snvId) => {
                  return this.snpDataStoreInstance.intToDnaSnv(snvId).snp_str;
                })
                .join(';') +
              '"'
          );
        } else if (key === 'gene_aa_snp_str') {
          vals.push(
            '"' +
              row[key]
                .map((snvId) => {
                  return this.snpDataStoreInstance.intToGeneAaSnv(snvId)
                    .snp_str;
                })
                .join(';') +
              '"'
          );
        } else if (key === 'protein_aa_snp_str') {
          vals.push(
            '"' +
              row[key]
                .map((snvId) => {
                  return this.snpDataStoreInstance.intToProteinAaSnv(snvId)
                    .snp_str;
                })
                .join(';') +
              '"'
          );
        } else {
          vals.push('"' + row[key] + '"');
        }
      });

      // Add location data
      vals = vals.concat(
        this.locationStoreInstance
          .getLocationStrFromId(row.location_id)
          .map((loc) => {
            return loc ? '"' + loc + '"' : loc;
          })
      );

      csvString.push(vals.join(','));
    });

    csvString = csvString.join('\n');

    const blob = new Blob([csvString]);
    const url = URL.createObjectURL(blob);
    downloadBlobURL(url, 'selected_sequence_metadata.csv');
  }

  @action
  downloadAggCaseData() {}

  @action
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

    downloadBlobURL(
      url,
      generateSelectionString(
        'data_agg_group_date',
        'csv',
        toJS(this.configStoreInstance.groupKey),
        toJS(this.configStoreInstance.dnaOrAa),
        toJS(this.configStoreInstance.selectedLocationNodes),
        toJS(this.configStoreInstance.dateRange)
      )
    );
  }

  @action
  downloadDataAggLocationGroupDate() {
    // console.log(this.dataAggLocationGroupDate);

    let locationData = JSON.parse(
      JSON.stringify(this.dataAggLocationGroupDate)
    );

    // Filter by date
    if (
      this.configStoreInstance.dateRange[0] != -1 ||
      this.configStoreInstance.dateRange[1] != -1
    ) {
      locationData = locationData.filter((row) => {
        return (
          (this.configStoreInstance.dateRange[0] == -1 ||
            row.date > this.configStoreInstance.dateRange[0]) &&
          (this.configStoreInstance.dateRange[1] == -1 ||
            row.date < this.configStoreInstance.dateRange[1])
        );
      });
    }

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

    downloadBlobURL(
      url,
      generateSelectionString(
        'data_agg_location_group_date',
        'csv',
        toJS(this.configStoreInstance.groupKey),
        toJS(this.configStoreInstance.dnaOrAa),
        toJS(this.configStoreInstance.selectedLocationNodes),
        toJS(this.configStoreInstance.dateRange)
      )
    );
  }

  @action
  downloadSnvFrequencies() {
    let csvString = 'snv,count\n';
    // this.groupCountArr.forEach((item) => {
    //   csvString += `${item[0]},${item[1]}\n`;
    // });

    const blob = new Blob([csvString]);
    const url = URL.createObjectURL(blob);

    downloadBlobURL(
      url,
      generateSelectionString(
        'snv_frequencies',
        'csv',
        toJS(this.configStoreInstance.groupKey),
        toJS(this.configStoreInstance.dnaOrAa),
        toJS(this.configStoreInstance.selectedLocationNodes),
        toJS(this.configStoreInstance.dateRange)
      )
    );
  }

  @action
  downloadSnvCooccurrence() {
    let csvString = 'combination,snv,count\n';
    this.snvCooccurrence.forEach((row) => {
      csvString += `${row.combi},${row.snv},${row.count}\n`;
    });

    const blob = new Blob([csvString]);
    const url = URL.createObjectURL(blob);

    downloadBlobURL(
      url,
      generateSelectionString(
        'snv_cooccurrence',
        'csv',
        toJS(this.configStoreInstance.groupKey),
        toJS(this.configStoreInstance.dnaOrAa),
        toJS(this.configStoreInstance.selectedLocationNodes),
        toJS(this.configStoreInstance.dateRange)
      )
    );
  }

  @action
  downloadCountryScoreData() {
    let jsonString = JSON.stringify(asyncDataStoreInstance.data.country_score);
    const blob = new Blob([jsonString]);
    const url = URL.createObjectURL(blob);

    downloadBlobURL(url, 'global_sequencing_coverage.json');
  }
}

export default DataStore;
