import { getLocationIds } from './location';
import { aggregate } from './transform';

import { DNA_OR_AA, COORDINATE_MODES, GROUPS } from '../constants/defs.json';
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

// dateRange is an array, [start, end]
function filterByDate(caseData, dateRange, dateKey = 'collection_date') {
  // Filter by date
  if (dateRange[0] > -1 && dateRange[1] > -1) {
    return caseData.filter((row) => {
      return row[dateKey] >= dateRange[0] && row[dateKey] <= dateRange[1];
    });
  }

  return caseData;
}

function getSnvFields({
  dnaOrAa,
  coordinateMode,
  selectedGroups,

  // SNV data
  intToDnaSnvMap,
  intToGeneAaSnvMap,
  intToProteinAaSnvMap,
  dnaSnvMap,
  geneAaSnvMap,
  proteinAaSnvMap,
}) {
  let selectedGroupIds;
  let snvEntry;
  let intToSnvMap;
  if (dnaOrAa === DNA_OR_AA.DNA) {
    selectedGroupIds = selectedGroups
      .map((item) => dnaSnvMap[item])
      .map((snpId) => (snpId === undefined ? -1 : snpId));
    snvEntry = 'dna_snp_str';
    intToSnvMap = intToDnaSnvMap;
  } else if (dnaOrAa === DNA_OR_AA.AA) {
    if (coordinateMode === COORDINATE_MODES.COORD_GENE) {
      selectedGroupIds = selectedGroups
        .map((item) => geneAaSnvMap[item])
        .map((snpId) => (snpId === undefined ? -1 : snpId));
      snvEntry = 'gene_aa_snp_str';
      intToSnvMap = intToGeneAaSnvMap;
    } else if (coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
      selectedGroupIds = selectedGroups
        .map((item) => proteinAaSnvMap[item])
        .map((snpId) => (snpId === undefined ? -1 : snpId));
      snvEntry = 'protein_aa_snp_str';
      intToSnvMap = intToProteinAaSnvMap;
    }
  }
  // Array to Set
  selectedGroupIds = new Set(selectedGroupIds);

  return {
    selectedGroupIds,
    snvEntry,
    intToSnvMap,
  };
}

function processSelectedSnvs({
  dnaOrAa,
  coordinateMode,
  selectedLocationNodes,
  countsPerLocation,
  filteredCaseData,
  selectedGroups,
  validGroups,

  // SNV data
  intToDnaSnvMap,
  intToGeneAaSnvMap,
  intToProteinAaSnvMap,
  dnaSnvMap,
  geneAaSnvMap,
  proteinAaSnvMap,
  snvColorMap,
}) {
  const { selectedGroupIds, snvEntry, intToSnvMap } = getSnvFields({
    dnaOrAa,
    coordinateMode,
    selectedGroups,

    // SNV data
    intToDnaSnvMap,
    intToGeneAaSnvMap,
    intToProteinAaSnvMap,
    dnaSnvMap,
    geneAaSnvMap,
    proteinAaSnvMap,
  });

  // If no SNVs are selected, then return empty arrays now
  // if (selectedGroupIds.size === 0) {
  //   return {
  //     dataAggLocationSnvDate: [],
  //     dataAggSnvDate: [],
  //     snvCooccurrence: [],
  //   };
  // }

  const dataAggLocationSnvDateObj = {};
  // Build a map of location_id --> node
  // Also while we're at it, create an entry for this node
  // in our data objects
  const selectedLocationIds = getLocationIds(selectedLocationNodes);
  const locationIdToNodeMap = {};
  for (let i = 0; i < selectedLocationNodes.length; i++) {
    selectedLocationIds[i].forEach((locationId) => {
      locationIdToNodeMap[locationId] = selectedLocationNodes[i].value;
    });
  }

  let _location;
  filteredCaseData.forEach((row) => {
    _location = locationIdToNodeMap[row.location_id];
    !(_location in dataAggLocationSnvDateObj) &&
      (dataAggLocationSnvDateObj[_location] = {});
    !(row.collection_date in dataAggLocationSnvDateObj[_location]) &&
      (dataAggLocationSnvDateObj[_location][row.collection_date] = {});

    // Check that every SNV ID is present
    // TODO: sort and then short-circuit check, that should be more efficient
    const matchingSnvIds = row[snvEntry].filter((snvId) =>
      selectedGroupIds.has(snvId)
    );
    let group;
    // "Other" group was selected, and this sequence has a SNV in "Other"
    if (
      selectedGroupIds.has(-1) &&
      row[snvEntry].some((id) => validGroups[id] === undefined)
    ) {
      group = GROUPS.OTHER_GROUP;
    } else if (
      matchingSnvIds.length != selectedGroupIds.size ||
      selectedGroups.length === 0
    ) {
      group = GROUPS.ALL_OTHER_GROUP;
    } else {
      group = matchingSnvIds.map((id) => intToSnvMap[id].snp_str).join(' + ');
    }

    !(group in dataAggLocationSnvDateObj[_location][row.collection_date]) &&
      (dataAggLocationSnvDateObj[_location][row.collection_date][group] = 0);
    dataAggLocationSnvDateObj[_location][row.collection_date][group] += 1;
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
          groupName: group
            .split(' + ')
            .map((group) => formatSnv(group, dnaOrAa))
            .join(' + '),
          cases_sum: dataAggLocationSnvDateObj[location][date][group],
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
  const dataAggSnvDate = aggregate({
    data: dataAggLocationSnvDate,
    groupby: ['date', 'group', 'groupName'],
    fields: ['cases_sum', 'color'],
    ops: ['sum', 'max'],
    as: ['cases_sum', 'color'],
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
  dnaOrAa,
  coordinateMode,
  selectedGroups,
  filteredCaseData,
  dateRange,

  // SNV data
  intToDnaSnvMap,
  intToGeneAaSnvMap,
  intToProteinAaSnvMap,
  dnaSnvMap,
  geneAaSnvMap,
  proteinAaSnvMap,
  snvColorMap,
}) {
  const { selectedGroupIds, snvEntry, intToSnvMap } = getSnvFields({
    dnaOrAa,
    coordinateMode,
    selectedGroups,

    // SNV data
    intToDnaSnvMap,
    intToGeneAaSnvMap,
    intToProteinAaSnvMap,
    dnaSnvMap,
    geneAaSnvMap,
    proteinAaSnvMap,
  });
  // Co-occurrence data

  // Re-filter by date
  filteredCaseData = filterByDate(filteredCaseData, dateRange);

  // Count SNVs for each single or combination of SNVs
  const snvCooccurrence = {};
  const selectedSnvCombinations = getCombinations(Array.from(selectedGroupIds));

  selectedSnvCombinations.forEach((combi) => {
    // Make an entry for this combination
    const combiKey = combi
      .map((snvId) => intToSnvMap[snvId].snp_str)
      .join(' + ');
    snvCooccurrence[combiKey] = {};

    // Loop thru all data, count SNVs
    filteredCaseData.forEach((row) => {
      if (!combi.every((snvId) => row[snvEntry].includes(snvId))) {
        return;
      }

      row[snvEntry].forEach((snvId) => {
        // Only check for SNVs that aren't in this combination
        if (!combi.includes(snvId)) {
          !(intToSnvMap[snvId].snp_str in snvCooccurrence[combiKey]) &&
            (snvCooccurrence[combiKey][intToSnvMap[snvId].snp_str] = 0);
          snvCooccurrence[combiKey][intToSnvMap[snvId].snp_str] += 1;
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
