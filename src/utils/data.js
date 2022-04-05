// Data processing functions

import { aggregate } from './transform';

import { LOW_FREQ_FILTER_TYPES, GROUP_MUTATION } from '../constants/defs.json';

/**
 * Collapse data aggregated by location, group, and date into
 * data by just group and date
 *
 * The input data will be an array of records, in this form
 * [{ location, date, group, count, etc ... }]
 * location = string of location name
 *
 * Resolve overlapping locations by identifying locations that are
 * strict subsets of other locations, and filtering their rows
 * out of the data, prior to aggregation
 */
export function removeSubsetLocations({
  selectedLocationNodes,
  aggLocationGroupDate,
}) {
  // We can use the location node "path", which is the
  // tree location embedded into a string, as a quick check for
  // strict subset-ness without having to iterate through and check
  // all child location IDs
  // First, sort nodes by their path's string length, which is a
  // quick approximation for the tree levels (nodes further down
  // the tree will have longer path string lengths)

  const compareByPathLength = (a, b) => {
    return a.path > b.path;
  };
  const sortedLocationNodes = selectedLocationNodes.sort(compareByPathLength);

  // Next, iterate through node list and prune nodes which are
  // subsets of existing nodes
  // This array is a list of location nodes that will be removed
  const keepLocationNodes = [];
  const prunedLocationNodes = [];
  sortedLocationNodes.forEach((node) => {
    // If this node is a subset of an existing pruned node, then skip
    // Otherwise, add it to the pruned list
    if (
      keepLocationNodes.some((existingNode) =>
        node.path.startsWith(existingNode.path)
      )
    ) {
      prunedLocationNodes.push(node);
    } else {
      keepLocationNodes.push(node);
    }
  });

  const prunedLocationNames = prunedLocationNodes.map((node) => node.value);

  // Filter out all pruned locations
  return aggLocationGroupDate.filter((record) => {
    return !prunedLocationNames.includes(record.location);
  });
}

/**
 * * The input data will be an array of records, in this form
 * [{ location, date, group, count, etc ... }]
 * location = string of location name
 *
 * And the goal is to collapse them into:
 * [{ date, group, count, etc ... }]
 *
 * In mutation mode, the `group_id` field will be an array of mutation ids
 * i.e., [1, 2, 3, ...]
 * In any other mode, the `group_id` field will be a string of the group name
 *
 * Aggregate the counts (summing them up) over the
 * date and group fields
 */
export function aggregateGroupDate({
  aggSequencesUniqueLocationGroupDate,
  groupKey,
}) {
  if (groupKey === GROUP_MUTATION) {
    // Same as below, except we first have to serialize our co-occurring mutations
    // into a hashable string, and then unpack it afterwards
    const aggGroupDate = aggregate({
      data: aggSequencesUniqueLocationGroupDate
        .filter((record) => {
          return record.group_id ? true : false;
        })
        .map((record) => {
          record.group_id_str = record.group_id.join(',');
          return record;
        }),
      groupby: ['collection_date', 'group_id_str'],
      fields: ['counts', 'group_id'],
      ops: ['sum', 'first'],
      as: ['counts', 'group_id'],
    });

    const aggGroup = aggregate({
      data: aggGroupDate,
      groupby: ['group_id_str'],
      fields: ['counts', 'group_id'],
      ops: ['sum', 'first'],
      as: ['counts', 'group_id'],
    });
    return {
      aggGroupDate,
      aggGroup,
    };
  } else {
    const aggGroupDate = aggregate({
      data: aggSequencesUniqueLocationGroupDate,
      groupby: ['collection_date', 'group_id'],
      fields: ['counts'],
      ops: ['sum'],
      as: ['counts'],
    });
    const aggGroup = aggregate({
      data: aggGroupDate,
      groupby: ['group_id'],
      fields: ['counts'],
      ops: ['sum'],
      as: ['counts'],
    });
    return {
      aggGroupDate,
      aggGroup,
    };
  }
}

/**
 * Given a list of location-date-group records
 * [{ location, date, group, count }, ...]
 * and low frequency collapse settings,
 * identify the groups to collapse into the "Other" group
 */
export function countGroups({ aggLocationGroupDate, groupKey }) {
  // Two main modes of data:
  // In mutation mode, the `group_id` is a list of mutation IDs
  // (co-occurring mutations, i.e., [1, 2, 3, ...])
  // In all other modes, `group_id` will be a string
  let aggSequencesGroup;
  if (groupKey === GROUP_MUTATION) {
    // Go through the list record by record and
    // tally up counts for each mutation ID
    const mutationIdCounts = new Map();
    aggLocationGroupDate.forEach((record) => {
      if (!record.group_id) return;
      record.group_id.forEach((mutationId) => {
        if (!mutationIdCounts.has(mutationId)) {
          mutationIdCounts.set(mutationId, 0);
        }
        mutationIdCounts.set(
          mutationId,
          mutationIdCounts.get(mutationId) + record.counts
        );
      });
    });

    // Convert Object to an array of records
    // [{ mutation_id, counts }]
    aggSequencesGroup = Array.from(mutationIdCounts.keys()).map(
      (mutationId) => {
        return {
          group_id: mutationId,
          counts: mutationIdCounts.get(mutationId),
        };
      }
    );
  } else {
    // Collapse data by group only, using the group_id string
    // as the groupby key
    aggSequencesGroup = aggregate({
      data: aggLocationGroupDate,
      groupby: ['group_id'],
      fields: ['counts'],
      ops: ['sum'],
      as: ['counts'],
    });
  }

  // Sort by counts in descending order
  return aggSequencesGroup.sort((a, b) => b.counts - a.counts);
}

/**
 * Use low frequency parameters to determine which groups to collapse
 */
export function getValidGroups({
  records,
  lowFreqFilterType,
  lowFreqFilterValue,
}) {
  let validGroups;
  if (lowFreqFilterType === LOW_FREQ_FILTER_TYPES.GROUP_COUNTS) {
    // Sort by counts in descending order
    validGroups = records
      .sort((a, b) => b.counts - a.counts)
      .slice(0, lowFreqFilterValue);
  } else if (lowFreqFilterType === LOW_FREQ_FILTER_TYPES.LOCAL_COUNTS) {
    validGroups = records.filter((group) => group.counts >= lowFreqFilterValue);
  }

  return validGroups.map((group) => group.group_id);
}

/**
 * Count the groups/mutations per location, and per location-date pair
 *
 * Input is structured as: [{ location, date, group, count }, ...]
 */
export function getLocationCounts({ aggLocationGroupDate }) {
  // Aggregate then calculate cumulative counts
  const countsPerLocationDate = aggregate({
    data: aggLocationGroupDate,
    groupby: ['location', 'collection_date'],
    fields: ['counts'],
    ops: ['sum'],
    as: ['counts'],
  }).sort((a, b) => {
    if (a.location === b.location) {
      return a.collection_date - b.collection_date;
    } else {
      return a.location > b.location;
    }
  });

  let cur_location =
    countsPerLocationDate.length === 0 ? '' : countsPerLocationDate[0].location;
  let cumulative_count = 0;
  countsPerLocationDate.forEach((record) => {
    if (cur_location !== record.location) {
      cumulative_count = 0;
      cur_location = record.location;
    }
    cumulative_count += record.counts;
    record.cumulative_count = cumulative_count;
  });

  // Convert into a map so we can do lookups like
  // countsPerLocationDate[location][date]
  const countsPerLocationDateMap = new Map();
  const cumulativeCountsPerLocationDateMap = new Map();
  countsPerLocationDate.forEach((record) => {
    if (!countsPerLocationDateMap.has(record.location)) {
      countsPerLocationDateMap.set(record.location, new Map());
      cumulativeCountsPerLocationDateMap.set(record.location, new Map());
    }
    countsPerLocationDateMap
      .get(record.location)
      .set(record.collection_date, record.counts);
    cumulativeCountsPerLocationDateMap
      .get(record.location)
      .set(record.collection_date, record.cumulative_count);
  });

  const countsPerLocationList = aggregate({
    data: countsPerLocationDate,
    groupby: ['location'],
    fields: ['counts'],
    ops: ['sum'],
    as: ['counts'],
  });
  // Convert array of records into key: val map
  const countsPerLocationMap = {};
  countsPerLocationList.forEach((record) => {
    countsPerLocationMap[record.location] = record.counts;
  });

  return {
    countsPerLocationDateMap,
    cumulativeCountsPerLocationDateMap,
    countsPerLocationMap,
  };
}

/**
 * Expand aggLocationGroupDate into
 * single mutation data
 * i.e., transform data where each row represents
 * a co-occurring mutation, into data where each row represents
 * individual mutations
 */
export function expandSingleMutationData({ aggLocationGroupDate }) {
  const singleMutationList = [];

  aggLocationGroupDate.forEach((record) => {
    if (!record.group_id) return;
    record.group_id.forEach((mutationId) => {
      singleMutationList.push({
        location: record.location,
        collection_date: record.collection_date,
        group_id: mutationId,
        counts: record.counts,
      });
    });
  });

  const aggLocationSingleMutationDate = aggregate({
    data: singleMutationList,
    groupby: ['location', 'collection_date', 'group_id'],
    fields: ['counts'],
    ops: ['sum'],
    as: ['counts'],
  });

  return aggLocationSingleMutationDate;
}
