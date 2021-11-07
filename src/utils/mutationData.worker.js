import { GROUPS } from '../constants/defs.json';
import { formatMutation } from './mutationUtils';

function getCombinations(arr) {
  var result = [];
  var f = function (prefix, arr) {
    for (var i = 0; i < arr.length; i++) {
      result.push([...prefix, arr[i]]);
      f([...prefix, arr[i]], arr.slice(i + 1));
    }
  };
  f([], arr);
  return result;
}

// mask unselected mutations
function processSelectedMutations({
  selectedGroupIds,
  intToMutationMap,
  dnaOrAa,
  countsPerLocationMap,
  // validGroups,
  aggLocationGroupDate,
  aggGroupDate,
  mutationColorMap,
}) {
  // If no mutations are selected, then return empty arrays now
  // if (selectedGroupIds.size === 0) {
  //   return {
  //     aggLocationSelectedMutationsDate: [],
  //     aggSelectedMutationsDate: [],
  //     mutationCooccurrence: [],
  //   };
  // }

  const aggLocationSelectedMutationsDateObj = {};

  const matchGroupName = Array.from(selectedGroupIds)
    .map((id) => intToMutationMap[id].snp_str)
    .join(' + ');

  // const validGroupMap = {};
  // validGroups.forEach((group) => {
  //   validGroupMap[group] = 1;
  // });

  aggLocationGroupDate.forEach((row) => {
    !(row.location in aggLocationSelectedMutationsDateObj) &&
      (aggLocationSelectedMutationsDateObj[row.location] = {});
    !(row.collection_date in aggLocationSelectedMutationsDateObj[row.location]) &&
      (aggLocationSelectedMutationsDateObj[row.location][row.collection_date] = {});

    // Check that every mutation ID is present
    // TODO: sort and then short-circuit check, that should be more efficient
    const matchingMutationIds = row.group_id.filter((mutationId) =>
      selectedGroupIds.has(mutationId)
    );
    let group;
    // "Other" group was selected, and this sequence has a mutation in "Other"
    if (
      selectedGroupIds.has(-1)
      // selectedGroupIds.has(-1) &&
      // row.group_id.some((id) => validGroupMap[id] === undefined)
    ) {
      group = GROUPS.OTHER_GROUP;
    } else if (
      matchingMutationIds.length != selectedGroupIds.size ||
      selectedGroupIds.size === 0
    ) {
      group = GROUPS.ALL_OTHER_GROUP;
    } else {
      group = matchGroupName;
    }

    !(
      group in aggLocationSelectedMutationsDateObj[row.location][row.collection_date]
    ) &&
      (aggLocationSelectedMutationsDateObj[row.location][row.collection_date][
        group
      ] = 0);
    aggLocationSelectedMutationsDateObj[row.location][row.collection_date][group] +=
      row.counts;
  });

  let groupKeys = [];
  const aggLocationSelectedMutationsDate = [];
  Object.keys(aggLocationSelectedMutationsDateObj).forEach((location) => {
    const dates = Object.keys(aggLocationSelectedMutationsDateObj[location]);
    dates.forEach((date) => {
      groupKeys = Object.keys(aggLocationSelectedMutationsDateObj[location][date]);
      groupKeys.forEach((group) => {
        aggLocationSelectedMutationsDate.push({
          location: location,
          collection_date: parseInt(date),
          group: group,
          group_name: group
            .split(' + ')
            .map((group) => formatMutation(group, dnaOrAa))
            .join(' + '),
          counts: aggLocationSelectedMutationsDateObj[location][date][group],
          // For multiple mutations, just use this nice blue color
          color:
            selectedGroupIds.size > 1 && group !== GROUPS.ALL_OTHER_GROUP
              ? '#07B'
              : mutationColorMap[group],
          location_counts: countsPerLocationMap[location],
        });
      });
    });
  });

  // Aggregate by mutation and date
  const aggSelectedMutationsDateObj = {};
  aggGroupDate.forEach((row) => {
    !(row.collection_date in aggSelectedMutationsDateObj) &&
      (aggSelectedMutationsDateObj[row.collection_date] = {});

    // Check that every mutation ID is present
    // TODO: sort and then short-circuit check, that should be more efficient
    const matchingMutationIds = row.group_id.filter((mutationId) =>
      selectedGroupIds.has(mutationId)
    );
    let group;
    // "Other" group was selected, and this sequence has a mutation in "Other"
    if (
      selectedGroupIds.has(-1) // &&
      // row.group_id.some((id) => validGroupMap[id] === undefined)
    ) {
      group = GROUPS.OTHER_GROUP;
    } else if (
      matchingMutationIds.length != selectedGroupIds.size ||
      selectedGroupIds.size === 0
    ) {
      group = GROUPS.ALL_OTHER_GROUP;
    } else {
      group = matchGroupName;
    }

    !(group in aggSelectedMutationsDateObj[row.collection_date]) &&
      (aggSelectedMutationsDateObj[row.collection_date][group] = 0);
    aggSelectedMutationsDateObj[row.collection_date][group] += row.counts;
  });

  const aggSelectedMutationsDate = [];
  Object.keys(aggSelectedMutationsDateObj).forEach((date) => {
    groupKeys = Object.keys(aggSelectedMutationsDateObj[date]);
    groupKeys.forEach((group) => {
      aggSelectedMutationsDate.push({
        collection_date: parseInt(date),
        group: group,
        group_name: group
          .split(' + ')
          .map((group) => formatMutation(group, dnaOrAa))
          .join(' + '),
        counts: aggSelectedMutationsDateObj[date][group],
        // For multiple mutations, just use this nice blue color
        color:
          selectedGroupIds.size > 1 && group !== GROUPS.ALL_OTHER_GROUP
            ? '#07B'
            : mutationColorMap[group],
      });
    });
  });
  // Convert dates back into ints
  aggSelectedMutationsDate.forEach((row) => {
    row.date = parseInt(row.date);
  });

  return {
    aggLocationSelectedMutationsDate,
    aggSelectedMutationsDate,
  };
}

function processCooccurrenceData({
  selectedGroupIds,
  intToMutationMap,
  dnaOrAa,
  aggSequencesGroup,
  // mutation data
  mutationColorMap,
}) {
  // Co-occurrence data

  // Count mutations for each single or combination of mutations
  const mutationCooccurrence = {};
  // Counts for each mutation combination - used to calculate fractions
  const mutationCombinationCounts = {};
  const selectedMutationCombinations = getCombinations(Array.from(selectedGroupIds));

  selectedMutationCombinations.forEach((combi) => {
    // Make an entry for this combination
    const combiKey = combi
      .map((mutationId) => intToMutationMap[mutationId].snp_str)
      .join(' + ');
    mutationCooccurrence[combiKey] = {};
    mutationCombinationCounts[combiKey] = 0;

    // Loop thru all data, count mutations
    aggSequencesGroup.forEach((row) => {
      if (!combi.every((mutationId) => row.group_id.includes(mutationId))) {
        return;
      }

      // Tally counts for this combination
      mutationCombinationCounts[combiKey] += row.counts;

      row.group_id.forEach((mutationId) => {
        // Only check for mutations that aren't in this combination
        if (!combi.includes(mutationId)) {
          !(intToMutationMap[mutationId].snp_str in mutationCooccurrence[combiKey]) &&
            (mutationCooccurrence[combiKey][intToMutationMap[mutationId].snp_str] = 0);
          mutationCooccurrence[combiKey][intToMutationMap[mutationId].snp_str] += row.counts;
        }
        // } else {
        //   !('None' in mutationCooccurrence[combiKey]) &&
        //     (mutationCooccurrence[combiKey]['None'] = 0);
        //   mutationCooccurrence[combiKey]['None'] += 1;
        // }
      });
    });
  });

  // Obj to list
  // This should also filter out any non-occurring combinations
  const mutationCooccurrenceList = [];
  Object.keys(mutationCooccurrence).forEach((combi) => {
    Object.keys(mutationCooccurrence[combi]).forEach((mutation) => {
      mutationCooccurrenceList.push({
        combi: combi,
        combiName: combi
          .split(' + ')
          .map((mut) => formatMutation(mut, dnaOrAa))
          .join(' + '),
        mutation: mut,
        mutationName: formatMutation(mutation, dnaOrAa),
        color: mutationColorMap[mutation],
        count: mutationCooccurrence[combi][mutation],
        fraction: mutationCooccurrence[combi][mutation] / mutationCombinationCounts[combi],
        combiCount: mutationCombinationCounts[combi],
      });
    });
  });
  // console.log(JSON.stringify(mutationCooccurrenceList));

  return {
    mutationCooccurrence: mutationCooccurrenceList,
  };
}

self.addEventListener(
  'message',
  function (e) {
    const data = e.data;
    //console.log('in downloadworker event listener', data);

    let result;
    if (data.type === 'processSelectedMutations') {
      result = processSelectedMutations(data);
    } else if (data.type === 'processCooccurrenceData') {
      result = processCooccurrenceData(data);
    }
    result.type = data.type;
    self.postMessage(result);
  },
  false
);
