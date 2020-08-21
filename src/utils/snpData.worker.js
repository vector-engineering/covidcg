import { DNA_OR_AA, COORDINATE_MODES } from '../constants/config';
import { getLocationIds } from './location';
import {
  dnaSnpToInt,
  geneAaSnpToInt,
  proteinAaSnpToInt,
  intToDnaSnp,
  intToGeneAaSnp,
  intToProteinAaSnp,
} from './snpData';

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

function processSelectedSnvs({
  dnaOrAa,
  coordinateMode,
  selectedLocationNodes,
  filteredCaseData,
  selectedGroups,
}) {
  let selectedGroupIds;
  let snvEntry;
  let intToSnvFunc;
  if (dnaOrAa === DNA_OR_AA.DNA) {
    selectedGroupIds = selectedGroups.map((item) => dnaSnpToInt(item));
    snvEntry = 'dna_snp_str';
    intToSnvFunc = intToDnaSnp;
  } else if (dnaOrAa === DNA_OR_AA.AA) {
    if (coordinateMode === COORDINATE_MODES.COORD_GENE) {
      selectedGroupIds = selectedGroups.map((item) => geneAaSnpToInt(item));
      snvEntry = 'gene_aa_snp_str';
      intToSnvFunc = intToGeneAaSnp;
    } else if (coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
      selectedGroupIds = selectedGroups.map((item) => proteinAaSnpToInt(item));
      snvEntry = 'protein_aa_snp_str';
      intToSnvFunc = intToProteinAaSnp;
    }
  }
  // Array to Set
  selectedGroupIds = new Set(selectedGroupIds);

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
    if (matchingSnvIds.length === 0) {
      group = 'Other';
    } else {
      group = matchingSnvIds.map((id) => intToSnvFunc(id).snp_str).join(' + ');
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
          cases_sum: dataAggLocationSnvDateObj[location][date][group],
          // color: getColorMethod(group),
        });
      });
    });
  });
  //console.log(dataAggLocationSnvDate);

  // Co-occurrence data
  // Count SNVs for each single or combination of SNVs
  const snvCooccurrence = {};
  const selectedSnvCombinations = getCombinations(Array.from(selectedGroupIds));

  selectedSnvCombinations.forEach((combi) => {
    // Make an entry for this combination
    const combiKey = combi
      .map((snvId) => intToSnvFunc(snvId).snp_str)
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
          !(intToSnvFunc(snvId).snp_str in snvCooccurrence[combiKey]) &&
            (snvCooccurrence[combiKey][intToSnvFunc(snvId).snp_str] = 0);
          snvCooccurrence[combiKey][intToSnvFunc(snvId).snp_str] += 1;
        } else {
          !('None' in snvCooccurrence[combiKey]) &&
            (snvCooccurrence[combiKey]['None'] = 0);
          snvCooccurrence[combiKey]['None'] += 1;
        }
      });
    });
  });
  //console.log(snvCooccurrence);

  return {
    dataAggLocationSnvDate,
    snvCooccurrence,
  };
}

self.addEventListener(
  'message',
  function (e) {
    const data = JSON.parse(e.data);
    //console.log('in downloadworker event listener', data);

    let result;
    if (data.type === 'processSelectedSnvs') {
      result = processSelectedSnvs(data);
    }
    // console.log(result);
    self.postMessage(JSON.stringify(result));
  },
  false
);
