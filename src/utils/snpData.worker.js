import { GROUPS } from '../constants/defs.json';
import { formatSnv } from './snpUtils';

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

// mask unselected SNVs
function processSelectedSnvs({
  selectedGroupIds,
  intToSnvMap,
  dnaOrAa,
  countsPerLocationMap,
  // validGroups,
  aggLocationGroupDate,
  aggGroupDate,
  snvColorMap,
}) {
  // If no SNVs are selected, then return empty arrays now
  // if (selectedGroupIds.size === 0) {
  //   return {
  //     aggLocationSelectedSnvsDate: [],
  //     aggSelectedSnvsDate: [],
  //     snvCooccurrence: [],
  //   };
  // }

  const aggLocationSelectedSnvsDateObj = {};

  const matchGroupName = Array.from(selectedGroupIds)
    .map((id) => intToSnvMap[id].snp_str)
    .join(' + ');

  // const validGroupMap = {};
  // validGroups.forEach((group) => {
  //   validGroupMap[group] = 1;
  // });

  aggLocationGroupDate.forEach((row) => {
    !(row.location in aggLocationSelectedSnvsDateObj) &&
      (aggLocationSelectedSnvsDateObj[row.location] = {});
    !(row.collection_date in aggLocationSelectedSnvsDateObj[row.location]) &&
      (aggLocationSelectedSnvsDateObj[row.location][row.collection_date] = {});

    // Check that every SNV ID is present
    // TODO: sort and then short-circuit check, that should be more efficient
    const matchingSnvIds = row.group_id.filter((snvId) =>
      selectedGroupIds.has(snvId)
    );
    let group;
    // "Other" group was selected, and this sequence has a SNV in "Other"
    if (
      selectedGroupIds.has(-1)
      // selectedGroupIds.has(-1) &&
      // row.group_id.some((id) => validGroupMap[id] === undefined)
    ) {
      group = GROUPS.OTHER_GROUP;
    } else if (
      matchingSnvIds.length != selectedGroupIds.size ||
      selectedGroupIds.size === 0
    ) {
      group = GROUPS.ALL_OTHER_GROUP;
    } else {
      group = matchGroupName;
    }

    !(
      group in aggLocationSelectedSnvsDateObj[row.location][row.collection_date]
    ) &&
      (aggLocationSelectedSnvsDateObj[row.location][row.collection_date][
        group
      ] = 0);
    aggLocationSelectedSnvsDateObj[row.location][row.collection_date][group] +=
      row.counts;
  });

  let groupKeys = [];
  const aggLocationSelectedSnvsDate = [];
  Object.keys(aggLocationSelectedSnvsDateObj).forEach((location) => {
    const dates = Object.keys(aggLocationSelectedSnvsDateObj[location]);
    dates.forEach((date) => {
      groupKeys = Object.keys(aggLocationSelectedSnvsDateObj[location][date]);
      groupKeys.forEach((group) => {
        aggLocationSelectedSnvsDate.push({
          location: location,
          collection_date: parseInt(date),
          group: group,
          group_name: group
            .split(' + ')
            .map((group) => formatSnv(group, dnaOrAa))
            .join(' + '),
          counts: aggLocationSelectedSnvsDateObj[location][date][group],
          // For multiple SNVs, just use this nice blue color
          color:
            selectedGroupIds.size > 1 && group !== GROUPS.ALL_OTHER_GROUP
              ? '#07B'
              : snvColorMap[group],
          location_counts: countsPerLocationMap[location],
        });
      });
    });
  });

  // Aggregate by SNV and date
  const aggSelectedSnvsDateObj = {};
  aggGroupDate.forEach((row) => {
    !(row.collection_date in aggSelectedSnvsDateObj) &&
      (aggSelectedSnvsDateObj[row.collection_date] = {});

    // Check that every SNV ID is present
    // TODO: sort and then short-circuit check, that should be more efficient
    const matchingSnvIds = row.group_id.filter((snvId) =>
      selectedGroupIds.has(snvId)
    );
    let group;
    // "Other" group was selected, and this sequence has a SNV in "Other"
    if (
      selectedGroupIds.has(-1) // &&
      // row.group_id.some((id) => validGroupMap[id] === undefined)
    ) {
      group = GROUPS.OTHER_GROUP;
    } else if (
      matchingSnvIds.length != selectedGroupIds.size ||
      selectedGroupIds.size === 0
    ) {
      group = GROUPS.ALL_OTHER_GROUP;
    } else {
      group = matchGroupName;
    }

    !(group in aggSelectedSnvsDateObj[row.collection_date]) &&
      (aggSelectedSnvsDateObj[row.collection_date][group] = 0);
    aggSelectedSnvsDateObj[row.collection_date][group] += row.counts;
  });

  const aggSelectedSnvsDate = [];
  Object.keys(aggSelectedSnvsDateObj).forEach((date) => {
    groupKeys = Object.keys(aggSelectedSnvsDateObj[date]);
    groupKeys.forEach((group) => {
      aggSelectedSnvsDate.push({
        collection_date: parseInt(date),
        group: group,
        group_name: group
          .split(' + ')
          .map((group) => formatSnv(group, dnaOrAa))
          .join(' + '),
        counts: aggSelectedSnvsDateObj[date][group],
        // For multiple SNVs, just use this nice blue color
        color:
          selectedGroupIds.size > 1 && group !== GROUPS.ALL_OTHER_GROUP
            ? '#07B'
            : snvColorMap[group],
      });
    });
  });
  // Convert dates back into ints
  aggSelectedSnvsDate.forEach((row) => {
    row.date = parseInt(row.date);
  });

  return {
    aggLocationSelectedSnvsDate,
    aggSelectedSnvsDate,
  };
}

function processCooccurrenceData({
  selectedGroupIds,
  intToSnvMap,
  dnaOrAa,
  aggSequencesGroup,
  // SNV data
  snvColorMap,
}) {
  // Co-occurrence data

  // Count SNVs for each single or combination of SNVs
  const snvCooccurrence = {};
  // Counts for each SNV combination - used to calculate fractions
  const snvCombinationCounts = {};
  const selectedSnvCombinations = getCombinations(Array.from(selectedGroupIds));

  selectedSnvCombinations.forEach((combi) => {
    // Make an entry for this combination
    const combiKey = combi
      .map((snvId) => intToSnvMap[snvId].snp_str)
      .join(' + ');
    snvCooccurrence[combiKey] = {};
    snvCombinationCounts[combiKey] = 0;

    // Loop thru all data, count SNVs
    aggSequencesGroup.forEach((row) => {
      if (!combi.every((snvId) => row.group_id.includes(snvId))) {
        return;
      }

      // Tally counts for this combination
      snvCombinationCounts[combiKey] += row.counts;

      row.group_id.forEach((snvId) => {
        // Only check for SNVs that aren't in this combination
        if (!combi.includes(snvId)) {
          !(intToSnvMap[snvId].snp_str in snvCooccurrence[combiKey]) &&
            (snvCooccurrence[combiKey][intToSnvMap[snvId].snp_str] = 0);
          snvCooccurrence[combiKey][intToSnvMap[snvId].snp_str] += row.counts;
        }
        // } else {
        //   !('None' in snvCooccurrence[combiKey]) &&
        //     (snvCooccurrence[combiKey]['None'] = 0);
        //   snvCooccurrence[combiKey]['None'] += 1;
        // }
      });
    });
  });

  // Obj to list
  // This should also filter out any non-occurring combinations
  const snvCooccurrenceList = [];
  Object.keys(snvCooccurrence).forEach((combi) => {
    Object.keys(snvCooccurrence[combi]).forEach((snv) => {
      snvCooccurrenceList.push({
        combi: combi,
        combiName: combi
          .split(' + ')
          .map((snv) => formatSnv(snv, dnaOrAa))
          .join(' + '),
        snv: snv,
        snvName: formatSnv(snv, dnaOrAa),
        color: snvColorMap[snv],
        count: snvCooccurrence[combi][snv],
        fraction: snvCooccurrence[combi][snv] / snvCombinationCounts[combi],
        combiCount: snvCombinationCounts[combi],
      });
    });
  });
  // console.log(JSON.stringify(snvCooccurrenceList));

  return {
    snvCooccurrence: snvCooccurrenceList,
  };
}

self.addEventListener(
  'message',
  function (e) {
    const data = e.data;
    //console.log('in downloadworker event listener', data);

    let result;
    if (data.type === 'processSelectedSnvs') {
      result = processSelectedSnvs(data);
    } else if (data.type === 'processCooccurrenceData') {
      result = processCooccurrenceData(data);
    }
    result.type = data.type;
    self.postMessage(result);
  },
  false
);
