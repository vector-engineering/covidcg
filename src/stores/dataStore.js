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
  downloadAcknowledgementsData,
  downloadAggCaseData,
} from '../utils/downloadWorkerWrapper';
import { aggregate } from '../utils/transform';
import { intToISO } from '../utils/date';
import { decryptAccessionIds } from '../utils/decrypt';
import { downloadBlobURL, generateSelectionString } from '../utils/download';
import { UIStoreInstance, configStoreInstance } from './rootStore';
import {
  LOW_FREQ_FILTER_TYPES,
  GROUP_KEYS,
  DNA_OR_AA,
  COORDINATE_MODES,
} from '../constants/config';
import { REFERENCE_GROUP } from '../constants/groups';
import { getGlobalGroupCounts } from '../utils/globalCounts';

import countryScoreData from '../../data/country_score.json';

const globalGroupCounts = getGlobalGroupCounts();

export const initialDataValues = {
  filteredCaseData: [],
  dataAggLocationGroupDate: [],
  dataAggGroupDate: [],
  dataAggGroup: [],
  changingPositions: {},
  groupsToKeep: [],
  groupCountArr: [],
  groupCountDateFilteredArr: [],

  // Metadata filtering
  numSequencesBeforeMetadataFiltering: 0,
  metadataCounts: {},
  metadataCountsAfterFiltering: {},

  selectedAccessionIds: [],
  selectedAckIds: [],

  dataAggLocationSnvDate: [],
  dataAggSnvDate: [],
  snvCooccurrence: [],

  countsPerLocation: {},
  countsPerLocationDate: {},
};

class ObservableDataStore {
  filteredCaseData = initialDataValues.filteredCaseData;
  dataAggLocationGroupDate = initialDataValues.dataAggLocationGroupDate;
  dataAggGroupDate = initialDataValues.dataAggGroupDate;
  dataAggGroup = initialDataValues.dataAggGroup;
  @observable changingPositions = initialDataValues.changingPositions;
  @observable groupsToKeep = initialDataValues.groupsToKeep;
  @observable groupCountArr = initialDataValues.groupCountArr;
  groupCountDateFilteredArr = initialDataValues.groupCountDateFilteredArr;
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

  countsPerLocation = initialDataValues.countsPerLocation;
  countsPerLocationDate = initialDataValues.countsPerLocationDate;

  constructor() {
    UIStoreInstance.onCaseDataStateStarted();
    this.updateCaseData();
  }

  @action
  updateGroupsToKeep() {
    if (this.groupCountArr === undefined) {
      return;
    }

    if (
      configStoreInstance.lowFreqFilterType ===
      LOW_FREQ_FILTER_TYPES.GROUP_COUNTS
    ) {
      this.groupsToKeep = this.groupCountArr
        .slice(0, configStoreInstance.maxGroupCounts)
        .map((item) => item[0]);
    } else if (
      configStoreInstance.lowFreqFilterType ===
      LOW_FREQ_FILTER_TYPES.LOCAL_COUNTS
    ) {
      this.groupsToKeep = this.groupCountArr
        .filter((item) => item[1] >= configStoreInstance.minLocalCountsToShow)
        .map((item) => item[0]);
    } else if (
      configStoreInstance.lowFreqFilterType ===
      LOW_FREQ_FILTER_TYPES.GLOBAL_COUNTS
    ) {
      let globalCounts;
      if (configStoreInstance.groupKey === GROUP_KEYS.GROUP_LINEAGE) {
        globalCounts = globalGroupCounts.lineage;
      } else if (configStoreInstance.groupKey === GROUP_KEYS.GROUP_CLADE) {
        globalCounts = globalGroupCounts.clade;
      } else if (configStoreInstance.groupKey === GROUP_KEYS.GROUP_SNV) {
        if (configStoreInstance.dnaOrAa === DNA_OR_AA.DNA) {
          globalCounts = globalGroupCounts.dna_snp;
        } else {
          if (
            configStoreInstance.coordinateMode === COORDINATE_MODES.COORD_GENE
          ) {
            globalCounts = globalGroupCounts.gene_aa_snp;
          } else if (
            configStoreInstance.coordinateMode ===
            COORDINATE_MODES.COORD_PROTEIN
          ) {
            globalCounts = globalGroupCounts.protein_aa_snp;
          }
        }
      }

      this.groupsToKeep = this.groupCountArr
        .filter(
          (item) =>
            globalCounts[item[0]] >= configStoreInstance.minGlobalCountsToShow
        )
        .map((item) => item[0]);
    }
    this.groupsToKeep.push(REFERENCE_GROUP);
  }

  @action
  updateAggCaseDataByGroup(callback) {
    UIStoreInstance.onAggCaseDataStarted();
    aggCaseDataByGroup(
      {
        filteredCaseData: JSON.parse(JSON.stringify(this.filteredCaseData)),
        dataAggGroupDate: JSON.parse(JSON.stringify(this.dataAggGroupDate)),
        coordinateMode: toJS(configStoreInstance.coordinateMode),
        coordinateRanges: toJS(configStoreInstance.getCoordinateRanges()),
        selectedGene: toJS(configStoreInstance.selectedGene),
        selectedProtein: toJS(configStoreInstance.selectedProtein),
        groupKey: toJS(configStoreInstance.groupKey),
        dnaOrAa: toJS(configStoreInstance.dnaOrAa),
        dateRange: toJS(configStoreInstance.dateRange),
      },
      ({
        dataAggGroup,
        changingPositions,
        groupCountArr,
        groupCountDateFilteredArr,
        selectedAccessionIds,
        selectedAckIds,
      }) => {
        // console.log(caseDataAggGroup);
        this.dataAggGroup = dataAggGroup;
        this.changingPositions = changingPositions;
        this.groupCountArr = groupCountArr;
        this.groupCountDateFilteredArr = groupCountDateFilteredArr;

        this.selectedAccessionIds = selectedAccessionIds;
        this.selectedAckIds = selectedAckIds;

        this.updateGroupsToKeep();
        // Fire callback
        if (callback !== undefined) {
          callback();
        }

        UIStoreInstance.onAggCaseDataFinished();
        UIStoreInstance.onCaseDataStateFinished();
      }
    );
  }

  @action
  emptyCaseData() {
    UIStoreInstance.onCaseDataStateStarted();
    UIStoreInstance.onAggCaseDataStarted();

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

    UIStoreInstance.onAggCaseDataFinished();
    UIStoreInstance.onCaseDataStateFinished();
  }

  @action
  updateCaseData(callback) {
    UIStoreInstance.onCaseDataStateStarted();
    processCaseData(
      {
        selectedLocationNodes: toJS(configStoreInstance.selectedLocationNodes),
        coordinateMode: toJS(configStoreInstance.coordinateMode),
        coordinateRanges: toJS(configStoreInstance.getCoordinateRanges()),
        selectedGene: toJS(configStoreInstance.selectedGene),
        selectedProtein: toJS(configStoreInstance.selectedProtein),
        groupKey: toJS(configStoreInstance.groupKey),
        dnaOrAa: toJS(configStoreInstance.dnaOrAa),
        selectedMetadataFields: toJS(
          configStoreInstance.selectedMetadataFields
        ),
        ageRange: toJS(configStoreInstance.ageRange),
        dateRange: toJS(configStoreInstance.dateRange),
        selectedGroups: toJS(configStoreInstance.selectedGroups),
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

        this.updateAggCaseDataByGroup(callback);
        if (configStoreInstance.groupKey === GROUP_KEYS.GROUP_SNV) {
          this.processSelectedSnvs();
        }
      }
    );
  }

  @action
  processSelectedSnvs() {
    UIStoreInstance.onSnvDataStarted();
    this.processCooccurrenceData();
    processSelectedSnvs(
      {
        dnaOrAa: toJS(configStoreInstance.dnaOrAa),
        coordinateMode: toJS(configStoreInstance.coordinateMode),
        selectedLocationNodes: toJS(configStoreInstance.selectedLocationNodes),
        filteredCaseData: this.filteredCaseData,
        selectedGroups: toJS(configStoreInstance.selectedGroups).map(
          (item) => item.group
        ),
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
    processCooccurrenceData(
      {
        dnaOrAa: toJS(configStoreInstance.dnaOrAa),
        coordinateMode: toJS(configStoreInstance.coordinateMode),
        filteredCaseData: this.filteredCaseData,
        selectedGroups: toJS(configStoreInstance.selectedGroups).map(
          (item) => item.group
        ),
        dateRange: toJS(configStoreInstance.dateRange),
      },
      ({ snvCooccurrence }) => {
        this.snvCooccurrence = snvCooccurrence;
        UIStoreInstance.onCooccurrenceDataFinished();
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
              toJS(configStoreInstance.groupKey),
              toJS(configStoreInstance.dnaOrAa),
              toJS(configStoreInstance.selectedLocationNodes),
              toJS(configStoreInstance.dateRange)
            )
          );
        }
      );
    });
  }

  @action
  downloadAcknowledgements() {
    // console.log('DOWNLOAD ACKNOWLEDGEMENTS');

    decryptAccessionIds(this.selectedAccessionIds).then((responseData) => {
      downloadAcknowledgementsData(
        {
          selectedAccessionIds: responseData['accession_ids'],
          selectedAckIds: this.selectedAckIds,
        },
        (res) => {
          // console.log(res);
          downloadBlobURL(
            res.blobURL,
            generateSelectionString(
              'acknowledgements',
              'csv',
              toJS(configStoreInstance.groupKey),
              toJS(configStoreInstance.dnaOrAa),
              toJS(configStoreInstance.selectedLocationNodes),
              toJS(configStoreInstance.dateRange)
            )
          );
        }
      );
    });
  }

  @action
  downloadAggCaseData() {
    downloadAggCaseData(
      {
        groupKey: configStoreInstance.groupKey,
        dnaOrAa: configStoreInstance.dnaOrAa,
        coordinateMode: configStoreInstance.coordinateMode,
        dataAggGroup: this.dataAggGroup,
      },
      (res) => {
        downloadBlobURL(
          res.blobURL,
          generateSelectionString(
            'agg_data',
            'csv',
            toJS(configStoreInstance.groupKey),
            toJS(configStoreInstance.dnaOrAa),
            toJS(configStoreInstance.selectedLocationNodes),
            toJS(configStoreInstance.dateRange)
          )
        );
      }
    );
  }

  @action
  downloadDataAggGroupDate() {
    // Write to a CSV string
    let csvString = `collection_date,${configStoreInstance.getGroupLabel()},count\n`;
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
        toJS(configStoreInstance.groupKey),
        toJS(configStoreInstance.dnaOrAa),
        toJS(configStoreInstance.selectedLocationNodes),
        toJS(configStoreInstance.dateRange)
      )
    );
  }

  @action
  downloadDataAggLocationGroupDate() {
    // console.log(this.dataAggLocationGroupDate);
    let csvString = `location,collection_date,${configStoreInstance.getGroupLabel()},count\n`;
    this.dataAggLocationGroupDate.forEach((row) => {
      csvString += `${row.location},${intToISO(parseInt(row.date))},${
        row.group
      },${row.cases_sum}\n`;
    });

    const blob = new Blob([csvString]);
    const url = URL.createObjectURL(blob);

    downloadBlobURL(
      url,
      generateSelectionString(
        'data_agg_location_group_date',
        'csv',
        toJS(configStoreInstance.groupKey),
        toJS(configStoreInstance.dnaOrAa),
        toJS(configStoreInstance.selectedLocationNodes),
        toJS(configStoreInstance.dateRange)
      )
    );
  }

  @action
  downloadSnvFrequencies() {
    let csvString = 'snv,count\n';
    this.groupCountArr.forEach((item) => {
      csvString += `${item[0]},${item[1]}\n`;
    });

    const blob = new Blob([csvString]);
    const url = URL.createObjectURL(blob);

    downloadBlobURL(
      url,
      generateSelectionString(
        'snv_frequencies',
        'csv',
        toJS(configStoreInstance.groupKey),
        toJS(configStoreInstance.dnaOrAa),
        toJS(configStoreInstance.selectedLocationNodes),
        toJS(configStoreInstance.dateRange)
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
        toJS(configStoreInstance.groupKey),
        toJS(configStoreInstance.dnaOrAa),
        toJS(configStoreInstance.selectedLocationNodes),
        toJS(configStoreInstance.dateRange)
      )
    );
  }

  @action
  downloadCountryScoreData() {
    let jsonString = JSON.stringify(countryScoreData);
    const blob = new Blob([jsonString]);
    const url = URL.createObjectURL(blob);

    downloadBlobURL(url, 'global_sequencing_coverage.json');
  }
}

export default ObservableDataStore;
