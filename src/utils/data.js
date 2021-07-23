// Data processing functions

import { aggregate } from './transform';
import { getLocationIdsByNode } from './location';

/**
 * Collapse data aggregated by location, group, and date into
 * data by just group and date
 *
 * The input data will be an array of records, in this form
 * [{ location, date, group, count, etc ... }]
 * location = string of location name
 *
 * And the goal is to collapse them into:
 * [{ date, group, count, etc ... }]
 *
 * By aggregating the counts (summing them up) over the
 * date and group fields
 *
 * Resolve overlapping locations by identifying locations that are
 * strict subsets of other locations, and filtering their rows
 * out of the data, prior to aggregation
 */
export function aggregateGroupDate({
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
  const prunedRecords = dataAggLocationGroupDate.filter((record) => {
    return !prunedLocationNames.includes(record.location);
  });

  return aggregate({
    data: prunedRecords,
    groupby: ['date', 'group'],
    fields: ['counts', 'color', 'group_name', 'group_id'],
    ops: ['sum', 'first', 'first', 'first'],
    as: ['counts', 'color', 'group_name', 'group_id'],
  });
}
