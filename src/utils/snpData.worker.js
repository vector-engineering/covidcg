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
  countsPerLocation,
  validGroups,
  aggSequencesLocationGroupDate,
  aggSequencesGroupDate,
  snvColorMap,
}) {
  // If no SNVs are selected, then return empty arrays now
  // if (selectedGroupIds.size === 0) {
  //   return {
  //     dataAggLocationSnvDate: [],
  //     dataAggSnvDate: [],
  //     snvCooccurrence: [],
  //   };
  // }

  const dataAggLocationSnvDateObj = {};

  const matchGroupName = Array.from(selectedGroupIds)
    .map((id) => intToSnvMap[id].snp_str)
    .join(' + ');

  const validGroupMap = {};
  validGroups.forEach((group) => {
    validGroupMap[group] = 1;
  });

  aggSequencesLocationGroupDate.forEach((row) => {
    !(row.location in dataAggLocationSnvDateObj) &&
      (dataAggLocationSnvDateObj[row.location] = {});
    !(row.collection_date in dataAggLocationSnvDateObj[row.location]) &&
      (dataAggLocationSnvDateObj[row.location][row.collection_date] = {});

    // Check that every SNV ID is present
    // TODO: sort and then short-circuit check, that should be more efficient
    const matchingSnvIds = row.group_id.filter((snvId) =>
      selectedGroupIds.has(snvId)
    );
    let group;
    // "Other" group was selected, and this sequence has a SNV in "Other"
    if (
      selectedGroupIds.has(-1) &&
      row.group_id.some((id) => validGroupMap[id] === undefined)
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

    !(group in dataAggLocationSnvDateObj[row.location][row.collection_date]) &&
      (dataAggLocationSnvDateObj[row.location][row.collection_date][group] = 0);
    dataAggLocationSnvDateObj[row.location][row.collection_date][group] +=
      row.counts;
  });

  let groupKeys = [];
  const dataAggLocationSnvDate = [];
  Object.keys(dataAggLocationSnvDateObj).forEach((location) => {
    const dates = Object.keys(dataAggLocationSnvDateObj[location]);
    dates.forEach((date) => {
      groupKeys = Object.keys(dataAggLocationSnvDateObj[location][date]);
      groupKeys.forEach((group) => {
        dataAggLocationSnvDate.push({
          location: location,
          date: parseInt(date),
          group: group,
          group_name: group
            .split(' + ')
            .map((group) => formatSnv(group, dnaOrAa))
            .join(' + '),
          counts: dataAggLocationSnvDateObj[location][date][group],
          // For multiple SNVs, just use this nice blue color
          color:
            selectedGroupIds.size > 1 && group !== GROUPS.ALL_OTHER_GROUP
              ? '#07B'
              : snvColorMap[group],
          location_counts: countsPerLocation[location],
        });
      });
    });
  });

  // Aggregate by SNV and date
  const dataAggSnvDateObj = {};
  aggSequencesGroupDate.forEach((row) => {
    !(row.collection_date in dataAggSnvDateObj) &&
      (dataAggSnvDateObj[row.collection_date] = {});

    // Check that every SNV ID is present
    // TODO: sort and then short-circuit check, that should be more efficient
    const matchingSnvIds = row.group_id.filter((snvId) =>
      selectedGroupIds.has(snvId)
    );
    let group;
    // "Other" group was selected, and this sequence has a SNV in "Other"
    if (
      selectedGroupIds.has(-1) &&
      row.group_id.some((id) => validGroupMap[id] === undefined)
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

    !(group in dataAggSnvDateObj[row.collection_date]) &&
      (dataAggSnvDateObj[row.collection_date][group] = 0);
    dataAggSnvDateObj[row.collection_date][group] += row.counts;
  });

  const dataAggSnvDate = [];
  Object.keys(dataAggSnvDateObj).forEach((date) => {
    groupKeys = Object.keys(dataAggSnvDateObj[date]);
    groupKeys.forEach((group) => {
      dataAggSnvDate.push({
        date: parseInt(date),
        group: group,
        group_name: group
          .split(' + ')
          .map((group) => formatSnv(group, dnaOrAa))
          .join(' + '),
        counts: dataAggSnvDateObj[date][group],
        // For multiple SNVs, just use this nice blue color
        color:
          selectedGroupIds.size > 1 && group !== GROUPS.ALL_OTHER_GROUP
            ? '#07B'
            : snvColorMap[group],
      });
    });
  });
  // Convert dates back into ints
  dataAggSnvDate.forEach((row) => {
    row.date = parseInt(row.date);
  });

  return {
    dataAggLocationSnvDate,
    dataAggSnvDate,
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
