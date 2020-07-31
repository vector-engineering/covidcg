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

class ObservableDataStore {
  @observable caseData = [];
  @observable changingPositions = {};
  @observable caseDataAggGroup = [];
  @observable selectedRows = [];
  @observable groupsToKeep = [];
  @observable lineageCountArr = {};

  // Metadata filtering
  @observable numSequencesBeforeMetadataFiltering = 0;
  @observable metadataCounts = {};

  // For location tab
  @observable aggLocationData = [];

  constructor() {
    UIStoreInstance.onCaseDataStateStarted();
    this.updateCaseData();
  }

  @action
  updateGroupsToKeep() {
    if (configStoreInstance.maxLineagesToShow > -1) {
      this.groupsToKeep = this.lineageCountArr
        .slice(0, configStoreInstance.maxLineagesToShow)
        .map((item) => item[0]);
    } else {
      this.groupsToKeep = this.lineageCountArr
        .filter((item) => item[1] >= configStoreInstance.minLocalCountsToShow)
        .map((item) => item[0]);
    }
    this.groupsToKeep.push('Reference');
  }

  @action
  updateAggCaseDataByGroup(suppressUIUpdate = false) {
    suppressUIUpdate ? null : UIStoreInstance.onAggCaseDataStarted();
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
      ({ caseDataAggGroup, changingPositions, lineageCountArr }) => {
        // console.log(caseDataAggGroup);
        this.caseDataAggGroup = caseDataAggGroup;
        this.changingPositions = changingPositions;
        this.lineageCountArr = lineageCountArr;
        // console.log('AGG_CASE_DATA FINISHED');

        this.updateGroupsToKeep();

        suppressUIUpdate ? null : UIStoreInstance.onAggCaseDataFinished();
        suppressUIUpdate ? null : UIStoreInstance.onCaseDataStateFinished();
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
  updateCaseData(suppressUIUpdate = false) {
    suppressUIUpdate ? null : UIStoreInstance.onCaseDataStateStarted();

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

        this.updateAggCaseDataByGroup((suppressUIUpdate = false));
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
