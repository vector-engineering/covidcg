// Data processing functions

import { aggregate } from './transform';

import { LOW_FREQ_FILTER_TYPES, GROUP_SNV } from '../constants/defs.json';

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
  dataAggLocationGroupDate,
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
  return dataAggLocationGroupDate.filter((record) => {
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
 * By aggregating the counts (summing them up) over the
 * date and group fields
 */
export function aggregateGroupDate({ dataAggLocationGroupDate }) {
  return aggregate({
    data: dataAggLocationGroupDate,
    groupby: ['date', 'group'],
    fields: ['counts', 'color', 'group_name', 'group_id'],
    ops: ['sum', 'first', 'first', 'first'],
    as: ['counts', 'color', 'group_name', 'group_id'],
  });
}

/**
 * Given a list of location-date-group records
 * [{ location, date, group, count }, ...]
 * and low frequency collapse settings,
 * identify the groups to collapse into the "Other" group
 */
export function countGroups({
  aggSequencesLocationGroupDate,
  groupKey,
  // lowFreqFilterType,
  // lowFreqFilterParams,
}) {
  // Two main modes of data:
  // In SNV mode, the `group_id` is a list of SNV IDs
  // (co-occurring mutations)
  // In all other modes, `group_id` will be a string
  let aggSequencesGroup;
  if (groupKey === GROUP_SNV) {
    // Go through the list record by record and
    // tally up counts for each SNV ID
    const snvIdCounts = {};
    aggSequencesLocationGroupDate.forEach((record) => {
      record.group_id.forEach((snvId) => {
        snvId = snvId.toString();
        if (!Object.prototype.hasOwnProperty.call(snvIdCounts, snvId)) {
          snvIdCounts[snvId] = 0;
        }
        snvIdCounts[snvId] += 1;
      });
    });

    // Convert Object to an array of records
    // [{ snv_id, counts }]
    aggSequencesGroup = Object.keys(snvIdCounts).map((snvIdStr) => {
      return { group_id: parseInt(snvIdStr), counts: snvIdCounts[snvIdStr] };
    });
  } else {
    // Collapse data by group only, using the group_id string
    // as the groupby key
    aggSequencesGroup = aggregate({
      data: aggSequencesLocationGroupDate,
      groupby: ['group_id'],
      fields: ['counts'],
      ops: ['sum'],
      as: ['counts'],
    });
  }

  // Sort by counts in descending order
  return aggSequencesGroup.sort((a, b) => b.counts - a.counts);

  // let validGroups;
  // if (lowFreqFilterType === LOW_FREQ_FILTER_TYPES.GROUP_COUNTS) {
  //   validGroups = aggSequencesGroup.slice(
  //     0,
  //     lowFreqFilterParams.maxGroupCounts
  //   );
  // } else if (lowFreqFilterType === LOW_FREQ_FILTER_TYPES.LOCAL_COUNTS) {
  //   validGroups = aggSequencesGroup.filter(
  //     (group) => group.count >= lowFreqFilterParams.minLocalCounts
  //   );
  // }

  // return {
  //   aggSequencesGroup,
  //   validGroups: validGroups.map((group) => group.group_id),
  // };
}
