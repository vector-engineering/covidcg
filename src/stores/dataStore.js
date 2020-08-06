import { observable, action, toJS } from 'mobx';
import {
  processCaseData,
  aggCaseDataByGroup,
} from '../utils/caseDataWorkerWrapper';
import {
  downloadAccessionIdsData,
  downloadAcknowledgementsData,
  downloadAggCaseData,
} from '../utils/downloadWorkerWrapper';
import { aggregate } from '../utils/transform';
import { decryptAccessionIds } from '../utils/decrypt';
import { downloadBlobURL, generateSelectionString } from '../utils/download';
import { UIStoreInstance, configStoreInstance } from './rootStore';
import {
  LOW_FREQ_FILTER_TYPES,
  GROUP_KEYS,
  DNA_OR_AA,
  COORDINATE_MODES,
} from '../constants/config';
import { getGlobalGroupCounts } from '../utils/globalCounts';

const globalGroupCounts = getGlobalGroupCounts();

export const initialDataValues = {
  dataAggLocationGroupDate: [],
  dataAggGroupDate: [],
  dataAggGroup: [],
  changingPositions: {},
  groupsToKeep: [],
  groupCountArr: [],

  // Metadata filtering
  numSequencesBeforeMetadataFiltering: 0,
  metadataCounts: {},

  selectedRowsHash: '',
  selectedRowsAndDateHash: '',
  selectedAccessionIds: [],
  selectedAckIds: [],
};

class ObservableDataStore {
  dataAggLocationGroupDate = initialDataValues.dataAggLocationGroupDate;
  dataAggGroupDate = initialDataValues.dataAggGroupDate;
  dataAggGroup = initialDataValues.dataAggGroup;
  @observable changingPositions = initialDataValues.changingPositions;
  @observable groupsToKeep = initialDataValues.groupsToKeep;
  @observable groupCountArr = initialDataValues.groupCountArr;
  @observable numSequencesBeforeMetadataFiltering =
    initialDataValues.numSequencesBeforeMetadataFiltering;
  @observable metadataCounts = initialDataValues.metadataCounts;

  @observable selectedRowsHash = initialDataValues.selectedRowsHash;
  @observable selectedRowsAndDateHash =
    initialDataValues.selectedRowsAndDateHash;
  selectedAccessionIds = initialDataValues.selectedAccessionIds;
  selectedAckIds = initialDataValues.selectedAckIds;

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
    this.groupsToKeep.push('Reference');
  }

  @action
  updateAggCaseDataByGroup(callback) {
    UIStoreInstance.onAggCaseDataStarted();
    aggCaseDataByGroup(
      {
        dataAggGroupDate: this.dataAggGroupDate,
        coordinateMode: toJS(configStoreInstance.coordinateMode),
        coordinateRanges: toJS(configStoreInstance.getCoordinateRanges()),
        selectedGene: toJS(configStoreInstance.selectedGene),
        selectedProtein: toJS(configStoreInstance.selectedProtein),
        groupKey: toJS(configStoreInstance.groupKey),
        dnaOrAa: toJS(configStoreInstance.dnaOrAa),
        dateRange: toJS(configStoreInstance.dateRange),
      },
      ({ dataAggGroup, changingPositions, groupCountArr }) => {
        // console.log(caseDataAggGroup);
        this.dataAggGroup = dataAggGroup;
        this.changingPositions = changingPositions;
        this.groupCountArr = groupCountArr;
        // console.log('AGG_CASE_DATA FINISHED');

        // Update hash for any listeners
        this.selectedRowsAndDateHash =
          this.selectedRowsHash +
          toJS(configStoreInstance.coordinateMode) +
          toJS(configStoreInstance.getCoordinateRanges()).toString() +
          configStoreInstance.selectedGene.gene +
          configStoreInstance.selectedProtein.protein +
          configStoreInstance.groupKey +
          configStoreInstance.dnaOrAa +
          toJS(configStoreInstance.dateRange)
            .map((date) => date.toString())
            .join(',');

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
    this.dataAggLocationGroupDate = [];
    this.dataAggGroupDate = [];
    this.dataAggGroup = [];
    this.changingPositions = {};
    this.selectedRowsHash = '';
    this.selectedAccessionIds = [];
    this.selectedAckIds = [];
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
      },
      ({
        aggCaseDataList,
        metadataCounts,
        numSequencesBeforeMetadataFiltering,
        selectedRowsHash,
        selectedAccessionIds,
        selectedAckIds,
      }) => {
        this.dataAggLocationGroupDate = aggCaseDataList;
        this.dataAggGroupDate = aggregate({
          data: aggCaseDataList,
          groupby: ['date', 'group'],
          fields: ['cases_sum', 'color'],
          ops: ['sum', 'max'],
          as: ['cases_sum', 'color'],
        });
        this.metadataCounts = metadataCounts;
        this.numSequencesBeforeMetadataFiltering = numSequencesBeforeMetadataFiltering;
        // console.log('CASE_DATA FINISHED');

        this.selectedRowsHash = selectedRowsHash;
        this.selectedAccessionIds = selectedAccessionIds;
        this.selectedAckIds = selectedAckIds;

        this.updateAggCaseDataByGroup(callback);
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
}

export default ObservableDataStore;
