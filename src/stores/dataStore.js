import { observable, action, toJS } from 'mobx';
import _ from 'underscore';

import {
  processCaseData,
  aggCaseDataByGroup,
} from '../utils/caseDataWorkerWrapper';
import {
  downloadAccessionIdsData,
  downloadAcknowledgementsData,
  downloadAggCaseData,
} from '../utils/downloadWorkerWrapper';
import { decryptAccessionIds } from '../utils/decrypt';
//import { getLineagesFromGene } from '../utils/lineageData';
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
  caseData: [],
  changingPositions: {},
  caseDataAggGroup: [],
  selectedRows: [],
  groupsToKeep: [],
  groupCountArr: {},

  // Metadata filtering
  numSequencesBeforeMetadataFiltering: 0,
  metadataCounts: {},

  // For location tab
  aggLocationData: [],
};

class ObservableDataStore {
  @observable caseData = initialDataValues.caseData;
  @observable changingPositions = initialDataValues.changingPositions;
  @observable caseDataAggGroup = initialDataValues.caseDataAggGroup;
  @observable selectedRows = initialDataValues.selectedRows;
  @observable groupsToKeep = initialDataValues.groupsToKeep;
  @observable groupCountArr = initialDataValues.groupCountArr;
  @observable numSequencesBeforeMetadataFiltering =
    initialDataValues.numSequencesBeforeMetadataFiltering;
  @observable metadataCounts = initialDataValues.metadataCounts;
  @observable aggLocationData = initialDataValues.aggLocationData;

  constructor() {
    UIStoreInstance.onCaseDataStateStarted();
    this.updateCaseData();
  }

  @action
  updateGroupsToKeep() {
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
        caseData: toJS(this.caseData),
        coordinateMode: toJS(configStoreInstance.coordinateMode),
        coordinateRanges: toJS(configStoreInstance.coordinateRanges),
        selectedGene: toJS(configStoreInstance.selectedGene),
        selectedProtein: toJS(configStoreInstance.selectedProtein),
        groupKey: toJS(configStoreInstance.groupKey),
        dnaOrAa: toJS(configStoreInstance.dnaOrAa),
        dateRange: toJS(configStoreInstance.dateRange),
      },
      ({ caseDataAggGroup, changingPositions, groupCountArr }) => {
        // console.log(caseDataAggGroup);
        this.caseDataAggGroup = caseDataAggGroup;
        this.changingPositions = changingPositions;
        this.groupCountArr = groupCountArr;
        // console.log('AGG_CASE_DATA FINISHED');

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
    this.caseData = [];
    this.selectedRows = [];
    this.caseDataAggGroup = [];
    this.changingPositions = {};
  }

  @action
  updateCaseData(callback) {
    UIStoreInstance.onCaseDataStateStarted();
    processCaseData(
      {
        selectedLocationIds: toJS(configStoreInstance.selectedLocationIds),
        coordinateMode: toJS(configStoreInstance.coordinateMode),
        coordinateRanges: toJS(configStoreInstance.coordinateRanges),
        selectedGene: toJS(configStoreInstance.selectedGene),
        selectedProtein: toJS(configStoreInstance.selectedProtein),
        groupKey: toJS(configStoreInstance.groupKey),
        dnaOrAa: toJS(configStoreInstance.dnaOrAa),
        selectedMetadataFields: toJS(
          configStoreInstance.selectedMetadataFields
        ),
        ageRange: toJS(configStoreInstance.ageRange),
        selectedLocationNodes: toJS(configStoreInstance.selectedLocationNodes),
      },
      ({
        aggCaseDataList,
        selectedRows,
        metadataCounts,
        numSequencesBeforeMetadataFiltering,
        aggLocationDataList,
      }) => {
        this.caseData = aggCaseDataList;
        this.selectedRows = selectedRows;
        this.metadataCounts = metadataCounts;
        this.numSequencesBeforeMetadataFiltering = numSequencesBeforeMetadataFiltering;
        this.aggLocationData = aggLocationDataList;
        // console.log('CASE_DATA FINISHED');

        this.updateAggCaseDataByGroup(callback);
      }
    );
  }

  @action
  downloadAccessionIds() {
    decryptAccessionIds(_.pluck(this.selectedRows, 'Accession ID')).then(
      (responseData) => {
        downloadAccessionIdsData(
          { accessionIds: responseData['accession_ids'] },
          (res) => {
            downloadBlobURL(
              res.blobURL,
              generateSelectionString(
                'accession_ids',
                'txt',
                configStoreInstance.groupKey,
                configStoreInstance.dnaOrAa,
                configStoreInstance.selectedLocationIds,
                configStoreInstance.dateRange
              )
            );
          }
        );
      }
    );
  }

  @action
  downloadAcknowledgements() {
    // console.log('DOWNLOAD ACKNOWLEDGEMENTS');

    decryptAccessionIds(_.pluck(this.selectedRows, 'Accession ID')).then(
      (responseData) => {
        // Make a deep copy of the selected rows
        let selectedRows = JSON.parse(JSON.stringify(toJS(this.selectedRows)));
        // Overwrite the existing hashed Accession IDs with the real ones
        for (let i = 0; i < selectedRows.length; i++) {
          selectedRows[i]['Accession ID'] = responseData['accession_ids'][i];
        }

        downloadAcknowledgementsData({ selectedRows }, (res) => {
          // console.log(res);
          downloadBlobURL(
            res.blobURL,
            generateSelectionString(
              'acknowledgements',
              'csv',
              configStoreInstance.groupKey,
              configStoreInstance.dnaOrAa,
              configStoreInstance.selectedLocationIds,
              configStoreInstance.dateRange
            )
          );
        });
      }
    );
  }

  @action
  downloadAggCaseData() {
    downloadAggCaseData(
      {
        groupKey: configStoreInstance.groupKey,
        dnaOrAa: configStoreInstance.dnaOrAa,
        coordinateMode: configStoreInstance.coordinateMode,
        caseDataAggGroup: toJS(this.caseDataAggGroup),
      },
      (res) => {
        downloadBlobURL(
          res.blobURL,
          generateSelectionString(
            'agg_data',
            'csv',
            configStoreInstance.groupKey,
            configStoreInstance.dnaOrAa,
            configStoreInstance.selectedLocationIds,
            configStoreInstance.dateRange
          )
        );
      }
    );
  }
}

export default ObservableDataStore;
