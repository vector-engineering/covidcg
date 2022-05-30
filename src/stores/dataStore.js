import { action, toJS } from 'mobx';
import { hostname } from '../config';
import {
  processSelectedMutations,
  processCooccurrenceData,
} from '../utils/mutationDataWorkerWrapper';
import { downloadBlobURL } from '../utils/download';
import { intToISO } from '../utils/date';
import { asyncDataStoreInstance } from '../components/App';
import { rootStoreInstance } from './rootStore';
import {
  GROUP_MUTATION,
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
  expandSingleMutationData,
} from '../utils/data';

export class DataStore {
  dataDate;
  numSequences;
  numSequencesAfterAllFiltering;
  timeToFetch;
  aggLocationGroupDate = [];
  aggGroupDate = [];
  aggLocationSingleMutationDate = [];

  aggLocationSelectedMutationsDate = [];
  aggSelectedMutationsDate = [];
  mutationCooccurrence = [];

  countsPerLocationDateMap = new Map();
  cumulativeCountsPerLocationDateMap = new Map();
  countsPerLocationMap = {};
  groupCounts = [];

  constructor() {}

  init() {
    this.dataDate = asyncDataStoreInstance.data.stats.data_date;
    this.numSequences = asyncDataStoreInstance.data.stats.num_sequences;

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

    const startTime = Date.now();

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
        ...rootStoreInstance.configStore.getSelectedLocations(),
        selected_reference: toJS(
          rootStoreInstance.configStore.selectedReference
        ),
        selected_metadata_fields:
          rootStoreInstance.configStore.getSelectedMetadataFields(),
        selected_group_fields: toJS(
          rootStoreInstance.configStore.selectedGroupFields
        ),
        ageRange: toJS(rootStoreInstance.configStore.ageRange),
        start_date: toJS(rootStoreInstance.configStore.startDate),
        end_date: toJS(rootStoreInstance.configStore.endDate),
        subm_start_date: toJS(rootStoreInstance.configStore.submStartDate),
        subm_end_date: toJS(rootStoreInstance.configStore.submEndDate),
      }),
    })
      .then((res) => {
        if (!res.ok) {
          throw res;
        }
        return res.json();
      })
      .then((res) => {
        // Time to fetch data, in seconds
        const endTime = Date.now();
        this.timeToFetch = (endTime - startTime) / 1000;

        this.aggLocationGroupDate = res['records'];
        // Convert time from unix epoch (seconds) to milliseconds
        this.aggLocationGroupDate.forEach((row) => {
          row.collection_date = row.collection_date * 1000;
        });

        // Pass coverage object onto mutationDataStore
        if (Object.prototype.hasOwnProperty.call(res, 'coverage')) {
          rootStoreInstance.mutationDataStore.coverage = res['coverage'];
        }

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

        // Aggregate, collapse locations
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
        // single mutation data
        // i.e., transform data where each row represents
        // a co-occurring mutation, into data where each row represents
        // individual mutations
        if (rootStoreInstance.configStore.groupKey === GROUP_MUTATION) {
          this.aggLocationSingleMutationDate = expandSingleMutationData({
            aggLocationGroupDate: this.aggLocationGroupDate,
          });
          // console.log(this.aggLocationSingleMutationDate);
        }

        rootStoreInstance.UIStore.onCaseDataStateFinished();

        if (rootStoreInstance.configStore.groupKey === GROUP_MUTATION) {
          this.processSelectedMutations();
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
  processSelectedMutations = () => {
    rootStoreInstance.UIStore.onMutationDataStarted();
    const { mutationColorMap } = rootStoreInstance.mutationDataStore;

    this.processCooccurrenceData();
    processSelectedMutations(
      {
        selectedGroupIds: rootStoreInstance.configStore.getSelectedGroupIds(),
        intToMutationMap: rootStoreInstance.configStore.getIntToMutationMap(),
        dnaOrAa: toJS(rootStoreInstance.configStore.dnaOrAa),
        countsPerLocationMap: this.countsPerLocationMap,
        // validGroups: this.validGroups,
        aggLocationGroupDate: this.aggLocationGroupDate,
        aggGroupDate: this.aggGroupDate,
        // Mutation data
        mutationColorMap,
      },
      ({ aggLocationSelectedMutationsDate, aggSelectedMutationsDate }) => {
        this.aggLocationSelectedMutationsDate =
          aggLocationSelectedMutationsDate;
        this.aggSelectedMutationsDate = aggSelectedMutationsDate;
        rootStoreInstance.UIStore.onMutationDataFinished();
      }
    );
  };

  @action
  processCooccurrenceData = () => {
    rootStoreInstance.UIStore.onCooccurrenceDataStarted();

    const { mutationColorMap } = rootStoreInstance.mutationDataStore;

    processCooccurrenceData(
      {
        selectedGroupIds: rootStoreInstance.configStore.getSelectedGroupIds(),
        intToMutationMap: rootStoreInstance.configStore.getIntToMutationMap(),
        dnaOrAa: toJS(rootStoreInstance.configStore.dnaOrAa),
        aggSequencesGroup: this.aggSequencesGroup,
        // mutation data
        mutationColorMap,
      },
      ({ mutationCooccurrence }) => {
        this.mutationCooccurrence = mutationCooccurrence;
        rootStoreInstance.UIStore.onCooccurrenceDataFinished();
      }
    );
  };

  @action
  downloadSelectedSequenceMetadata = ({
    selectedFields,
    mutationFormat,
    selectedReference,
  }) => {
    rootStoreInstance.UIStore.onDownloadStarted();

    fetch(hostname + '/download_metadata', {
      method: 'POST',
      headers: {
        Accept: 'text/csv',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        ...rootStoreInstance.configStore.getSelectedLocations(),
        selected_reference: selectedReference,
        selected_group_fields: toJS(
          rootStoreInstance.configStore.selectedGroupFields
        ),
        selected_metadata_fields:
          rootStoreInstance.configStore.getSelectedMetadataFields(),
        ageRange: toJS(rootStoreInstance.configStore.ageRange),
        start_date: toJS(rootStoreInstance.configStore.startDate),
        end_date: toJS(rootStoreInstance.configStore.endDate),
        subm_start_date: toJS(rootStoreInstance.configStore.submStartDate),
        subm_end_date: toJS(rootStoreInstance.configStore.submEndDate),
        // Pass an array of only the fields that were selected
        selected_fields: Object.keys(selectedFields).filter(
          (field) => selectedFields[field]
        ),
        mutation_format: mutationFormat,
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

  // TODO:
  // We should probably change this request to use form data
  // or query params, so that we can open the request in a new
  // window. Right now the gzipped file has to be loaded as a
  // blob on the front-end, and this will get unsustainable
  // if the user tries to download 100,000+ genomes
  @action
  downloadGenomes = ({ compress }) => {
    rootStoreInstance.UIStore.onDownloadStarted();

    fetch(hostname + '/download_genomes', {
      method: 'POST',
      headers: {
        Accept: 'application/gzip',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        ...rootStoreInstance.configStore.getSelectedLocations(),
        selected_group_fields: toJS(
          rootStoreInstance.configStore.selectedGroupFields
        ),
        selected_metadata_fields:
          rootStoreInstance.configStore.getSelectedMetadataFields(),
        ageRange: toJS(rootStoreInstance.configStore.ageRange),
        start_date: toJS(rootStoreInstance.configStore.startDate),
        end_date: toJS(rootStoreInstance.configStore.endDate),
        subm_start_date: toJS(rootStoreInstance.configStore.submStartDate),
        subm_end_date: toJS(rootStoreInstance.configStore.submEndDate),
        compress,
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
        downloadBlobURL(url, compress ? 'genomes.fa.gz' : 'genomes.fa');
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

    let intToMutationMap = rootStoreInstance.configStore.getIntToMutationMap();
    this.aggLocationGroupDate.forEach((row) => {
      // Location data, date
      csvString += `"${row.location}","${intToISO(row.collection_date)}",`;

      // Group
      if (rootStoreInstance.configStore.groupKey === GROUP_MUTATION) {
        csvString += `"${row.group_id
          .map((mutationId) => intToMutationMap[mutationId].mutation_str)
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

  downloadAggGroup() {
    let csvString = `${rootStoreInstance.configStore.getGroupLabel()},count,percent\n`;
    let intToMutationMap = rootStoreInstance.configStore.getIntToMutationMap();

    this.groupCounts.forEach((row) => {
      // Group
      if (rootStoreInstance.configStore.groupKey === GROUP_MUTATION) {
        csvString += `"${intToMutationMap[row.group_id].mutation_str}",`;
      } else {
        csvString += `"${row.group_id}",`;
      }

      // Counts
      csvString += `${row.counts},`;
      // Percent
      csvString += `${row.counts / this.numSequencesAfterAllFiltering}\n`;
    });

    const blob = new Blob([csvString]);
    const url = URL.createObjectURL(blob);

    downloadBlobURL(url, `group_counts.csv`);
  }

  downloadAggGroupDate() {
    // Write to a CSV string
    let csvString = `collection_date,"${rootStoreInstance.configStore.getGroupLabel()}",count\n`;
    this.aggGroupDate.forEach((row) => {
      let groupName;
      // If we're in mutation mode, then we have to map mutation IDs back
      // to a co-occurring mutation string representation
      if (rootStoreInstance.configStore.groupKey === GROUP_MUTATION) {
        groupName = row.group_id
          .map(
            (mutationId) =>
              rootStoreInstance.mutationDataStore.intToMutation(
                rootStoreInstance.configStore.dnaOrAa,
                rootStoreInstance.configStore.coordinateMode,
                mutationId
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
    if (rootStoreInstance.configStore.groupKey === GROUP_MUTATION) {
      locationData = toJS(this.aggLocationSingleMutationDate);
      // Get mutation data
      locationData.forEach((record) => {
        let mutation = rootStoreInstance.mutationDataStore.intToMutation(
          rootStoreInstance.configStore.dnaOrAa,
          rootStoreInstance.configStore.coordinateMode,
          record.group_id
        );
        record.group_name = mutation.name;
        record.group = mutation.mutation_str;
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

  downloadMutationFrequencies() {
    let csvString = 'mutation,';
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

    let mutation;

    // groupCounts is a list of records
    // [{group_id: int, group: str, counts: int, color: str, group_name: str}]
    this.groupCounts
      .sort((a, b) => b.counts - a.counts) // Sort by counts, descending order
      .forEach((record) => {
        mutation = rootStoreInstance.mutationDataStore.intToMutation(
          rootStoreInstance.configStore.dnaOrAa,
          rootStoreInstance.configStore.coordinateMode,
          record.group_id
        );

        // Add mutation fields
        if (mutation.mutation_str === GROUPS.REFERENCE_GROUP) {
          csvString += mutation.mutation_str + ',';
          csvString += fields.slice().fill('').join(',');
        } else {
          csvString +=
            mutation.name +
            ',' +
            fields.map((field) => mutation[field]).join(',');
        }
        // Add counts
        csvString += `,${record.counts}\n`;
      });
    // console.log(csvString);

    const blob = new Blob([csvString]);
    const url = URL.createObjectURL(blob);

    downloadBlobURL(url, 'mutation_frequencies.csv');
  }

  downloadMutationCooccurrence() {
    let csvString = 'selected_mutations,cooccurs_with,';
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

    // Have to convert mutation string into integer, then into mutation object
    const mutationToIntMap =
      rootStoreInstance.configStore.getMutationToIntMap();
    const intToMutationMap =
      rootStoreInstance.configStore.getIntToMutationMap();
    let mutation;

    this.mutationCooccurrence
      .sort((a, b) => {
        // Group by the mutation combination, then sort by counts, descending order
        if (a.combi < b.combi) {
          return -1;
        } else if (a.combi > b.combi) {
          return 1;
        } else {
          return b.count - a.count;
        }
      })
      .forEach((row) => {
        // Add combination and co-occuring mutation name
        csvString += `${row.combiName},${row.mutationName},`;
        // Add mutation fields
        mutation = intToMutationMap[mutationToIntMap[row.mutation]];
        csvString += fields.map((field) => mutation[field]).join(',');

        // Add counts and fraction
        csvString += `,${row.count},${row.combiCount},${row.fraction}\n`;
      });

    const blob = new Blob([csvString]);
    const url = URL.createObjectURL(blob);

    downloadBlobURL(url, 'mutation_cooccurrence.csv');
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
