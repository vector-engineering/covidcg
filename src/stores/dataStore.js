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
import { GROUP_KEYS } from '../constants/config';
import { asyncDataStoreInstance } from '../components/App';
import { rootStoreInstance } from './rootStore';

//import countryScoreData from '../../data/country_score.json';

const { UIStore, configStore } = rootStoreInstance || {};

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

  countsPerLocation: {},
  countsPerLocationDate: {},
  validGroups: {},
  countsPerGroup: {},
  countsPerGroupDateFiltered: [],
};

class ObservableDataStore {
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

  countsPerLocation = initialDataValues.countsPerLocation;
  countsPerLocationDate = initialDataValues.countsPerLocationDate;
  validGroups = initialDataValues.validGroups;
  countsPerGroup = initialDataValues.countsPerGroup;
  countsPerGroupDateFiltered = initialDataValues.countsPerGroupDateFiltered;

  constructor() {
    UIStore.onCaseDataStateStarted();
  }

  @action
  updateAggCaseDataByGroup(callback) {
    UIStore.onAggCaseDataStarted();
    aggCaseDataByGroup(
      {
        filteredCaseData: JSON.parse(JSON.stringify(this.filteredCaseData)),
        dataAggGroupDate: JSON.parse(JSON.stringify(this.dataAggGroupDate)),
        coordinateMode: toJS(configStore.coordinateMode),
        coordinateRanges: toJS(configStore.getCoordinateRanges()),
        selectedGene: toJS(configStore.selectedGene),
        selectedProtein: toJS(configStore.selectedProtein),
        groupKey: toJS(configStore.groupKey),
        dnaOrAa: toJS(configStore.dnaOrAa),
        dateRange: toJS(configStore.dateRange),
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

        UIStore.onAggCaseDataFinished();
        UIStore.onCaseDataStateFinished();
      }
    );
  }

  @action
  emptyCaseData() {
    UIStore.onCaseDataStateStarted();
    UIStore.onAggCaseDataStarted();

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

    UIStore.onAggCaseDataFinished();
    UIStore.onCaseDataStateFinished();
  }

  @action
  updateCaseData(callback) {
    UIStore.onCaseDataStateStarted();
    processCaseData(
      {
        selectedLocationNodes: toJS(configStore.selectedLocationNodes),
        groupKey: toJS(configStore.groupKey),
        dnaOrAa: toJS(configStore.dnaOrAa),
        selectedMetadataFields: toJS(configStore.selectedMetadataFields),
        ageRange: toJS(configStore.ageRange),
        dateRange: toJS(configStore.dateRange),
        selectedGroups: toJS(configStore.selectedGroups),

        coordinateMode: toJS(configStore.coordinateMode),
        coordinateRanges: toJS(configStore.getCoordinateRanges()),
        selectedGene: toJS(configStore.selectedGene),
        selectedProtein: toJS(configStore.selectedProtein),

        lowFreqFilterType: configStore.lowFreqFilterType,
        maxGroupCounts: configStore.maxGroupCounts,
        minLocalCountsToShow: configStore.minLocalCountsToShow,
        minGlobalCountsToShow: configStore.minGlobalCountsToShow,
        globalGroupCounts: toJS(asyncDataStoreInstance.globalGroupCounts),
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
        if (configStore.groupKey === GROUP_KEYS.GROUP_SNV) {
          this.processSelectedSnvs();
        }
      }
    );
  }

  @action
  processSelectedSnvs() {
    UIStore.onSnvDataStarted();
    this.processCooccurrenceData();
    processSelectedSnvs(
      {
        dnaOrAa: toJS(configStore.dnaOrAa),
        coordinateMode: toJS(configStore.coordinateMode),
        selectedLocationNodes: toJS(configStore.selectedLocationNodes),
        countsPerLocation: this.countsPerLocation,
        filteredCaseData: this.filteredCaseData,
        selectedGroups: toJS(configStore.selectedGroups).map(
          (item) => item.group
        ),
        validGroups: this.validGroups,
      },
      ({ dataAggLocationSnvDate, dataAggSnvDate }) => {
        this.dataAggLocationSnvDate = dataAggLocationSnvDate;
        this.dataAggSnvDate = dataAggSnvDate;
        UIStore.onSnvDataFinished();
      }
    );
  }

  @action
  processCooccurrenceData() {
    UIStore.onCooccurrenceDataStarted();
    processCooccurrenceData(
      {
        dnaOrAa: toJS(configStore.dnaOrAa),
        coordinateMode: toJS(configStore.coordinateMode),
        filteredCaseData: this.filteredCaseData,
        selectedGroups: toJS(configStore.selectedGroups).map(
          (item) => item.group
        ),
        dateRange: toJS(configStore.dateRange),
      },
      ({ snvCooccurrence }) => {
        this.snvCooccurrence = snvCooccurrence;
        UIStore.onCooccurrenceDataFinished();
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
              toJS(configStore.groupKey),
              toJS(configStore.dnaOrAa),
              toJS(configStore.selectedLocationNodes),
              toJS(configStore.dateRange)
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
              toJS(configStore.groupKey),
              toJS(configStore.dnaOrAa),
              toJS(configStore.selectedLocationNodes),
              toJS(configStore.dateRange)
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
        groupKey: configStore.groupKey,
        dnaOrAa: configStore.dnaOrAa,
        coordinateMode: configStore.coordinateMode,
        dataAggGroup: this.dataAggGroup,
      },
      (res) => {
        downloadBlobURL(
          res.blobURL,
          generateSelectionString(
            'agg_data',
            'csv',
            toJS(configStore.groupKey),
            toJS(configStore.dnaOrAa),
            toJS(configStore.selectedLocationNodes),
            toJS(configStore.dateRange)
          )
        );
      }
    );
  }

  @action
  downloadDataAggGroupDate() {
    // Write to a CSV string
    let csvString = `collection_date,${configStore.getGroupLabel()},count\n`;
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
        toJS(configStore.groupKey),
        toJS(configStore.dnaOrAa),
        toJS(configStore.selectedLocationNodes),
        toJS(configStore.dateRange)
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
    if (configStore.dateRange[0] != -1 || configStore.dateRange[1] != -1) {
      locationData = locationData.filter((row) => {
        return (
          (configStore.dateRange[0] == -1 ||
            row.date > configStore.dateRange[0]) &&
          (configStore.dateRange[1] == -1 ||
            row.date < configStore.dateRange[1])
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

    let csvString = `location,collection_date,${configStore.getGroupLabel()},count,location_date_count\n`;

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
        toJS(configStore.groupKey),
        toJS(configStore.dnaOrAa),
        toJS(configStore.selectedLocationNodes),
        toJS(configStore.dateRange)
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
        toJS(configStore.groupKey),
        toJS(configStore.dnaOrAa),
        toJS(configStore.selectedLocationNodes),
        toJS(configStore.dateRange)
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
        toJS(configStore.groupKey),
        toJS(configStore.dnaOrAa),
        toJS(configStore.selectedLocationNodes),
        toJS(configStore.dateRange)
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

export default ObservableDataStore;
