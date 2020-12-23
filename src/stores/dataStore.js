import { observable, action, toJS } from 'mobx';
import {
  processCaseData,
  aggCaseDataByGroup,
} from '../utils/caseDataWorkerWrapper';
import {
  processSelectedSnvs,
  processCooccurrenceData,
} from '../utils/snpDataWorkerWrapper';
import {
  downloadAccessionIdsData,
  // downloadAcknowledgementsData,
  downloadAggCaseData,
} from '../utils/downloadWorkerWrapper';
import { aggregate } from '../utils/transform';
import { intToISO } from '../utils/date';
import { decryptAccessionIds } from '../utils/decrypt';
import { downloadBlobURL, generateSelectionString } from '../utils/download';
import { GROUP_SNV } from '../constants/config';
import { asyncDataStoreInstance } from '../components/App';
import { rootStoreInstance } from './rootStore';

//import countryScoreData from '../../data/country_score.json';

export const initialDataValues = {
  filteredCaseData: [],
  dataAggLocationGroupDate: [],
  dataAggGroupDate: [],
  dataAggGroup: [],
  changingPositions: {},

  // Metadata filtering
  numSequencesBeforeMetadataFiltering: 0,
  metadataCounts: {},
  metadataCountsAfterFiltering: {},

  selectedAccessionIds: [],
  selectedAckIds: [],

  dataAggLocationSnvDate: [],
  dataAggSnvDate: [],
  snvCooccurrence: [],

  globalGroupCounts: {},
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
  groupDataStoreInstance;

  dataDate;
  numSequences;
  rawCaseData = [];
  filteredCaseData = initialDataValues.filteredCaseData;
  dataAggLocationGroupDate = initialDataValues.dataAggLocationGroupDate;
  dataAggGroupDate = initialDataValues.dataAggGroupDate;
  dataAggGroup = initialDataValues.dataAggGroup;
  @observable changingPositions = initialDataValues.changingPositions;
  @observable numSequencesBeforeMetadataFiltering =
    initialDataValues.numSequencesBeforeMetadataFiltering;
  @observable metadataCounts = initialDataValues.metadataCounts;
  @observable metadataCountsAfterFiltering =
    initialDataValues.metadataCountsAfterFiltering;

  selectedAccessionIds = initialDataValues.selectedAccessionIds;
  selectedAckIds = initialDataValues.selectedAckIds;

  dataAggLocationSnvDate = initialDataValues.dataAggLocationSnvDate;
  dataAggSnvDate = initialDataValues.dataAggSnvDate;
  snvCooccurrence = initialDataValues.snvCooccurrence;

  globalGroupCounts = initialDataValues.globalGroupCounts;
  countsPerLocation = initialDataValues.countsPerLocation;
  countsPerLocationDate = initialDataValues.countsPerLocationDate;
  validGroups = initialDataValues.validGroups;
  countsPerGroup = initialDataValues.countsPerGroup;
  countsPerGroupDateFiltered = initialDataValues.countsPerGroupDateFiltered;

  constructor() {}

  init() {
    this.dataDate = asyncDataStoreInstance.data.data_date;
    this.rawCaseData = asyncDataStoreInstance.data.case_data.map((row) => {
      row.collection_date = new Date(row.collection_date).getTime();
      return row;
    });
    this.rawCaseData = this.rawCaseData.filter((row) => {
      // Remove cases before 2019-12-15 and after the dataDate
      return !(
        row.collection_date < 1576368000000 ||
        row.collection_date > new Date(this.dataDate).getTime()
      );
    });
    this.numSequences = this.rawCaseData.length;

    // Calculate global group counts
    // Make a copy
    const globalGroupCounts = Object.assign(
      {},
      asyncDataStoreInstance.data.global_group_counts
    );

    // Replace integer IDs with SNP strings
    Object.keys(globalGroupCounts.dna_snp).forEach((snpId) => {
      globalGroupCounts.dna_snp[snpId.toString()] =
        globalGroupCounts.dna_snp[snpId];
    });
    Object.keys(globalGroupCounts.gene_aa_snp).forEach((snpId) => {
      globalGroupCounts.gene_aa_snp[snpId.toString()] =
        globalGroupCounts.gene_aa_snp[snpId];
    });
    Object.keys(globalGroupCounts.protein_aa_snp).forEach((snpId) => {
      globalGroupCounts.protein_aa_snp[snpId.toString()] =
        globalGroupCounts.protein_aa_snp[snpId];
    });
    this.globalGroupCounts = globalGroupCounts;

    this.UIStoreInstance = rootStoreInstance.UIStore;
    this.configStoreInstance = rootStoreInstance.configStore;
    this.snpDataStoreInstance = rootStoreInstance.snpDataStore;
    this.groupDataStoreInstance = rootStoreInstance.groupDataStore;

    this.updateCaseData();
  }

  @action
  updateAggCaseDataByGroup(callback) {
    this.UIStoreInstance.onAggCaseDataStarted();

    const {
      snvColorMap,
      intToDnaSnvMap,
      intToGeneAaSnvMap,
      intToProteinAaSnvMap,
    } = this.snpDataStoreInstance;

    const { groupSnvMap, groupColorMap } = this.groupDataStoreInstance;

    aggCaseDataByGroup(
      {
        filteredCaseData: JSON.parse(JSON.stringify(this.filteredCaseData)),
        dataAggGroupDate: JSON.parse(JSON.stringify(this.dataAggGroupDate)),
        coordinateMode: toJS(this.configStoreInstance.coordinateMode),
        coordinateRanges: toJS(this.configStoreInstance.getCoordinateRanges()),
        selectedGene: toJS(this.configStoreInstance.selectedGene),
        selectedProtein: toJS(this.configStoreInstance.selectedProtein),
        groupKey: toJS(this.configStoreInstance.groupKey),
        dnaOrAa: toJS(this.configStoreInstance.dnaOrAa),
        dateRange: toJS(this.configStoreInstance.dateRange),

        // SNV data
        snvColorMap,
        intToDnaSnvMap,
        intToGeneAaSnvMap,
        intToProteinAaSnvMap,

        // Lineage data
        groupSnvMap,
        groupColorMap,
      },
      ({
        dataAggGroup,
        changingPositions,
        selectedAccessionIds,
        selectedAckIds,
        countsPerGroupDateFiltered,
      }) => {
        // console.log(caseDataAggGroup);
        this.dataAggGroup = dataAggGroup;
        this.changingPositions = changingPositions;
        this.countsPerGroupDateFiltered = countsPerGroupDateFiltered;

        this.selectedAccessionIds = selectedAccessionIds;
        this.selectedAckIds = selectedAckIds;

        // Fire callback
        if (callback !== undefined) {
          callback();
        }

        this.UIStoreInstance.onAggCaseDataFinished();
        this.UIStoreInstance.onCaseDataStateFinished();
      }
    );
  }

  @action
  emptyCaseData() {
    this.UIStoreInstance.onCaseDataStateStarted();
    this.UIStoreInstance.onAggCaseDataStarted();

    this.dataAggLocationGroupDate = initialDataValues.dataAggLocationGroupDate;
    this.dataAggGroupDate = initialDataValues.dataAggGroupDate;
    this.dataAggGroup = initialDataValues.dataAggGroup;
    this.changingPositions = initialDataValues.changingPositions;

    this.selectedAccessionIds = initialDataValues.selectedAccessionIds;
    this.selectedAckIds = initialDataValues.selectedAckIds;

    this.dataAggLocationSnvDate = initialDataValues.dataAggLocationSnvDate;
    this.dataAggSnvDate = initialDataValues.dataAggSnvDate;
    this.snvCooccurrence = initialDataValues.snvCooccurrence;

    this.countsPerLocation = initialDataValues.countsPerLocation;
    this.countsPerLocationDate = initialDataValues.countsPerLocationDate;

    this.UIStoreInstance.onAggCaseDataFinished();
    this.UIStoreInstance.onCaseDataStateFinished();
  }

  @action
  updateCaseData(callback) {
    this.UIStoreInstance.onCaseDataStateStarted();

    const {
      intToDnaSnvMap,
      intToGeneAaSnvMap,
      intToProteinAaSnvMap,
      snvColorMap,
    } = this.snpDataStoreInstance;
    const { groupSnvMap, groupColorMap } = this.groupDataStoreInstance;

    processCaseData(
      {
        rawCaseData: this.rawCaseData,

        selectedLocationNodes: toJS(
          this.configStoreInstance.selectedLocationNodes
        ),
        groupKey: toJS(this.configStoreInstance.groupKey),
        dnaOrAa: toJS(this.configStoreInstance.dnaOrAa),
        selectedMetadataFields: toJS(
          this.configStoreInstance.selectedMetadataFields
        ),
        ageRange: toJS(this.configStoreInstance.ageRange),
        dateRange: toJS(this.configStoreInstance.dateRange),
        selectedGroups: toJS(this.configStoreInstance.selectedGroups),

        coordinateMode: toJS(this.configStoreInstance.coordinateMode),
        coordinateRanges: toJS(this.configStoreInstance.getCoordinateRanges()),
        selectedGene: toJS(this.configStoreInstance.selectedGene),
        selectedProtein: toJS(this.configStoreInstance.selectedProtein),

        lowFreqFilterType: this.configStoreInstance.lowFreqFilterType,
        maxGroupCounts: this.configStoreInstance.maxGroupCounts,
        minLocalCountsToShow: this.configStoreInstance.minLocalCountsToShow,
        minGlobalCountsToShow: this.configStoreInstance.minGlobalCountsToShow,
        globalGroupCounts: this.globalGroupCounts,

        // SNV data
        intToDnaSnvMap,
        intToGeneAaSnvMap,
        intToProteinAaSnvMap,
        snvColorMap,

        // Lineage data
        groupSnvMap,
        groupColorMap,
      },
      ({
        filteredCaseData,
        dataAggLocationGroupDate,
        dataAggGroupDate,
        metadataCounts,
        metadataCountsAfterFiltering,
        numSequencesBeforeMetadataFiltering,
        countsPerLocation,
        countsPerLocationDate,
        validGroups,
        countsPerGroup,
      }) => {
        this.filteredCaseData = filteredCaseData;
        this.dataAggLocationGroupDate = dataAggLocationGroupDate;
        this.dataAggGroupDate = dataAggGroupDate;
        this.metadataCounts = metadataCounts;
        this.metadataCountsAfterFiltering = metadataCountsAfterFiltering;
        this.numSequencesBeforeMetadataFiltering = numSequencesBeforeMetadataFiltering;
        // console.log('CASE_DATA FINISHED');

        this.countsPerLocation = countsPerLocation;
        this.countsPerLocationDate = countsPerLocationDate;
        this.validGroups = validGroups;
        this.countsPerGroup = countsPerGroup;

        this.updateAggCaseDataByGroup(callback);
        if (this.configStoreInstance.groupKey === GROUP_SNV) {
          this.processSelectedSnvs();
        }
      }
    );
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
  downloadAccessionIds() {
    decryptAccessionIds(this.selectedAccessionIds).then((responseData) => {
      downloadAccessionIdsData(
        { accessionIds: responseData['accession_ids'] },
        (res) => {
          downloadBlobURL(
            res.blobURL,
            generateSelectionString(
              'accession_ids',
              'txt',
              toJS(this.configStoreInstance.groupKey),
              toJS(this.configStoreInstance.dnaOrAa),
              toJS(this.configStoreInstance.selectedLocationNodes),
              toJS(this.configStoreInstance.dateRange)
            )
          );
        }
      );
    });
  }

  // @action
  // downloadAcknowledgements() {
  //   // console.log('DOWNLOAD ACKNOWLEDGEMENTS');

  //   decryptAccessionIds(this.selectedAccessionIds).then((responseData) => {
  //     downloadAcknowledgementsData(
  //       {
  //         selectedAccessionIds: responseData['accession_ids'],
  //         selectedAckIds: this.selectedAckIds,
  //       },
  //       (res) => {
  //         // console.log(res);
  //         downloadBlobURL(
  //           res.blobURL,
  //           generateSelectionString(
  //             'acknowledgements',
  //             'csv',
  //             toJS(configStoreInstance.groupKey),
  //             toJS(configStoreInstance.dnaOrAa),
  //             toJS(configStoreInstance.selectedLocationNodes),
  //             toJS(configStoreInstance.dateRange)
  //           )
  //         );
  //       }
  //     );
  //   });
  // }

  @action
  downloadAggCaseData() {
    const {
      intToDnaSnvMap,
      intToGeneAaSnvMap,
      intToProteinAaSnvMap,
    } = this.snpDataStoreInstance;

    const { groupSnvMap, groupColorMap } = this.groupDataStoreInstance;

    downloadAggCaseData(
      {
        groupKey: this.configStoreInstance.groupKey,
        dnaOrAa: this.configStoreInstance.dnaOrAa,
        coordinateMode: this.configStoreInstance.coordinateMode,
        dataAggGroup: this.dataAggGroup,

        // SNV data
        intToDnaSnvMap,
        intToGeneAaSnvMap,
        intToProteinAaSnvMap,

        // Lineage data
        groupSnvMap,
        groupColorMap,
      },
      (res) => {
        downloadBlobURL(
          res.blobURL,
          generateSelectionString(
            'agg_data',
            'csv',
            toJS(this.configStoreInstance.groupKey),
            toJS(this.configStoreInstance.dnaOrAa),
            toJS(this.configStoreInstance.selectedLocationNodes),
            toJS(this.configStoreInstance.dateRange)
          )
        );
      }
    );
  }

  @action
  downloadDataAggGroupDate() {
    // Write to a CSV string
    let csvString = `collection_date,${this.configStoreInstance.getGroupLabel()},count\n`;
    this.dataAggGroupDate.forEach((row) => {
      csvString += `${intToISO(parseInt(row.date))},${row.group},${
        row.cases_sum
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
      groupby: ['location', 'date', 'group', 'groupName'],
      fields: ['cases_sum', 'location_counts'],
      ops: ['sum', 'max'],
      as: ['cases_sum', 'location_counts'],
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
      },${row.cases_sum},${row.location_date_count}\n`;
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
    let jsonString = JSON.stringify(
      asyncDataStoreInstance.data.countryScoreData
    );
    const blob = new Blob([jsonString]);
    const url = URL.createObjectURL(blob);

    downloadBlobURL(url, 'global_sequencing_coverage.json');
  }
}

export default DataStore;
