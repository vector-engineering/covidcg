import initialCaseData from '../../data/case_data2.json';
import { intToDnaSnp, intToGeneAaSnp, intToProteinAaSnp } from './snpData';
import {
  getDnaSnpsFromGroup,
  getGeneAaSnpsFromGroup,
  getProteinAaSnpsFromGroup,
} from './lineageData';
import {
  warmColors,
  coolColors,
  snpColorArray,
  cladeColorArray,
  incrementColor,
} from '../constants/colors';
import { countMetadataFields } from './metadata';
import _ from 'underscore';

const processedCaseData = _.reject(
  _.map(initialCaseData, (row) => {
    row.collection_date = new Date(row.collection_date).getTime();
    return row;
  }),
  (row) => {
    // Remove cases before 2019-12-15
    return row.collection_date < 1576368000000;
  }
);

let _warmColors = Object.assign({}, warmColors);
let _coolColors = Object.assign({}, coolColors);
let _snpColorArray = [...snpColorArray];
let _cladeColorArray = [...cladeColorArray];

function iterateAndGetColor(colorobj, iter, i = 0) {
  if (i < iter) {
    return iterateAndGetColor(colorobj.child, iter, i + 1);
  } else {
    colorobj.colors.push(incrementColor(colorobj.colors.shift(), 1));
    return colorobj.colors[0];
  }
}

const getLineageColor = _.memoize((group) => {
  const match = group.match(/\./g);
  const dots = match ? match.length : 0;
  if (group.charAt(0) === 'A') {
    return iterateAndGetColor(_warmColors, dots);
  } else if (group.charAt(0) === 'B') {
    return iterateAndGetColor(_coolColors, dots);
  }
});

const getSnpColor = _.memoize(() => {
  _snpColorArray.push(incrementColor(_snpColorArray.shift(), 1));
  return _snpColorArray[0];
});

const getCladeColor = _.memoize(() => {
  _cladeColorArray.push(incrementColor(_cladeColorArray.shift(), 1));
  return _cladeColorArray[0];
});

function convertToObj(list) {
  const obj = {};
  list.forEach((item) => {
    obj[item] = 1;
  });
  return obj;
}

function filterByLocation(caseData, locationIds) {
  // Skip if no locations selected
  if (locationIds.length === 0) {
    return caseData;
  }

  const locationObj = convertToObj(
    // Location IDs are a list of lists (one list for each selected node),
    // so collapse into just one list
    locationIds.reduce((memo, list) => memo.concat(list), [])
  );

  return caseData.filter((row) => {
    return locationObj[row.location_id] === 1;
  });
}

function getParent(group) {
  let result = group.split('.');
  if (result.length < 1) {
    return 'root_node';
  }
  result.pop();

  result = result.join('.');
  return result;
}

function filterByCoordinateRange({ caseData, coordinateRanges }) {
  let newCaseData = [];

  caseData.forEach((row) => {
    // Only keep SNPs that are within
    row['dna_snp_str'] = _.filter(row['dna_snp_str'], (snpId) => {
      let snpObj = intToDnaSnp(snpId);
      // Keep the SNP if it falls into any one of the ranges
      return _.some(coordinateRanges, (range) => {
        return snpObj.pos >= range[0] && snpObj.pos <= range[1];
      });
    });
    newCaseData.push(row);
  });

  return newCaseData;
}

function filterByGeneOrProtein({
  caseData,
  coordinateMode,
  selectedGene,
  selectedProtein,
}) {
  let newCaseData = [];

  if (coordinateMode === 'gene') {
    // Don't need to do this if we selected all genes
    if (selectedGene.gene === 'All Genes') {
      return caseData;
    }

    caseData.forEach((row) => {
      row['gene_aa_snp_str'] = _.filter(row['gene_aa_snp_str'], (snpId) => {
        let snpObj = intToGeneAaSnp(snpId);
        return snpObj.gene === selectedGene.gene;
      });
      newCaseData.push(row);
    });
  } else if (coordinateMode === 'protein') {
    // Don't need to do this if we selected all proteins
    if (selectedProtein.protein === 'All Proteins') {
      return caseData;
    }

    caseData.forEach((row) => {
      row['protein_aa_snp_str'] = _.filter(
        row['protein_aa_snp_str'],
        (snpId) => {
          let snpObj = intToProteinAaSnp(snpId);
          return snpObj.protein === selectedProtein.protein;
        }
      );
      newCaseData.push(row);
    });
  } else {
    return caseData;
  }

  return newCaseData;
}

// dateRange is an array, [start, end]
function filterByDate(caseData, dateRange) {
  // Filter by date
  if (dateRange[0] > -1 && dateRange[1] > -1) {
    return _.filter(caseData, (row) => {
      return row.date >= dateRange[0] && row.date <= dateRange[1];
    });
  }

  return caseData;
}

function filterByMetadataFieldsAndAgeRange(
  caseData,
  selectedMetadataFields,
  ageRange
) {
  // Remove any fields from the selectedMetadataFields object that are empty arrays
  let metadataFields = {};
  Object.keys(selectedMetadataFields).forEach((field) => {
    if (selectedMetadataFields[field].length > 0) {
      metadataFields[field] = selectedMetadataFields[field];
    }
  });

  let remove;
  caseData = _.reject(caseData, (row) => {
    remove = false;
    Object.keys(metadataFields).forEach((field) => {
      if (!metadataFields[field].includes(row[field])) {
        remove = true;
      }
    });

    // Age range
    if (
      (ageRange[0] !== null && row['age_start'] < ageRange[0]) ||
      (ageRange[1] !== null && row['age_end'] > ageRange[1])
    ) {
      remove = true;
    }

    return remove;
  });

  return caseData;
}

function getGroupKeys(row, groupKey, dnaOrAa, coordinateMode) {
  // If we're grouping by lineage or SNP signature, then the group
  // keys are literal
  if (groupKey === 'lineage') {
    return [row['lineage']];
  } else if (groupKey === 'clade') {
    return [row['clade']];
  }
  // If we're grouping by SNP, then the value is a semicolon-delimited
  // list of SNPs that we should treat separately
  else if (groupKey === 'snp') {
    if (dnaOrAa === 'dna') {
      return row['dna_snp_str'];
    } else {
      if (coordinateMode === 'gene') {
        return row['gene_aa_snp_str'];
      } else if (coordinateMode === 'protein') {
        return row['protein_aa_snp_str'];
      }
    }
  }
}

function processCaseData({
  selectedLocationIds,
  coordinateMode,
  coordinateRanges,
  selectedGene,
  selectedProtein,
  groupKey,
  dnaOrAa,
  selectedMetadataFields,
  ageRange,
  selectedGroups,
  selectedLocationNodes,
}) {
  // let caseData = _.map(_caseData, (row) => Object.assign({}, row));
  let caseData = JSON.parse(JSON.stringify(processedCaseData));
  //console.log('filtering by selectedLocationIds', selectedLocationIds);

  // Filter by location
  caseData = filterByLocation(caseData, selectedLocationIds);
  // console.log(caseData.length, 'rows remaining after location filtering');
  // Filter by coordinate range (DNA mode only)
  if (groupKey === 'snp' && dnaOrAa === 'dna') {
    caseData = filterByCoordinateRange({ caseData, coordinateRanges });
    // console.log(
    //   caseData.length,
    //   'rows remaining after coordinate range filtering'
    // );
  }
  // Filter by gene/protein (AA mode only, gene or protein selected)
  if (groupKey === 'snp' && dnaOrAa === 'aa') {
    caseData = filterByGeneOrProtein({
      caseData,
      coordinateMode,
      selectedGene,
      selectedProtein,
    });
    // console.log(caseData.length, 'rows remaining after gene/protein filtering');
  }

  // Get the initial number of sequences, prior to metadata filtering
  const numSequencesBeforeMetadataFiltering = caseData.length;

  const metadataCounts = countMetadataFields(caseData);

  // Filter by metadata fields and age range
  caseData = filterByMetadataFieldsAndAgeRange(
    caseData,
    selectedMetadataFields,
    ageRange
  );
  // console.log(caseData.length, 'rows remaining after metadata filtering');

  // Get a list of Accession IDs and sample dates that are currently selected
  let selectedRows = _.map(caseData, (row) => {
    return row;
  });

  // For the location tab:
  const aggLocationData = {};
  const aggLocationDataByDate = {}; // Group by location and date
  const aggLocationDataByGroup = {}; // Group by location and group (lineage/clade/snp)
  const countsPerLocation = {};

  // Build a map of location_id --> node
  // Also while we're at it, create an entry for this node
  // in our data objects
  const locationIdToNodeMap = {};
  for (let i = 0; i < selectedLocationNodes.length; i++) {
    selectedLocationIds[i].forEach((locationId) => {
      locationIdToNodeMap[locationId] = selectedLocationNodes[i].value;
    });
    aggLocationDataByDate[selectedLocationNodes[i].value] = {};
    aggLocationDataByGroup[selectedLocationNodes[i].value] = {};
    countsPerLocation[selectedLocationNodes[i].value] = 0;
  }

  // If we have selected groups (lineages/snps/clades), then filter for that
  // const selectedGroupKeys = _.pluck(selectedGroups, 'group');
  let groupKeys = [];
  let location;
  caseData.forEach((row) => {
    countsPerLocation[locationIdToNodeMap[row.location_id]] += 1;
    groupKeys = getGroupKeys(row, groupKey, dnaOrAa, coordinateMode);

    location = locationIdToNodeMap[row.location_id];
    !(location in aggLocationData) && (aggLocationData[location] = {});
    !(row.collection_date in aggLocationData[location]) &&
      (aggLocationData[location][row.collection_date] = {});
    groupKeys.forEach((group) => {
      !(group in aggLocationData[location][row.collection_date]) &&
        (aggLocationData[location][row.collection_date][group] = 0);
      aggLocationData[location][row.collection_date][group] += 1;
    });
  });

  // Compute total counts per group
  // Object.keys(aggLocationDataByDate).forEach((location) => {
  //   countsPerLocation[location] = Object.values(
  //     aggLocationDataByDate[location]
  //   ).reduce((memo, val) => memo + val, 0);
  // });

  let getColorMethod;
  if (groupKey === 'lineage') {
    getColorMethod = getLineageColor;
  } else if (groupKey === 'clade') {
    getColorMethod = getCladeColor;
  } else if (groupKey === 'snp') {
    getColorMethod = getSnpColor;
  }

  const aggLocationDataList = [];
  Object.keys(aggLocationData).forEach((location) => {
    const dates = Object.keys(aggLocationData[location]);
    dates.forEach((date) => {
      groupKeys = Object.keys(aggLocationData[location][date]);
      groupKeys.forEach((group) => {
        aggLocationDataList.push({
          location: location,
          date: parseInt(date),
          group: group,
          cases_sum: aggLocationData[location][date][group],
          location_counts: countsPerLocation[location],
          color: getColorMethod(group),
        });
      });
    });
  });
  // console.log(JSON.stringify(aggLocationDataList));

  // Group by grouping key and sample date
  // console.log('Grouping by', groupKey, 'and sample date');
  const aggCaseData = {};
  let row = {};
  // let groupKeys = [];
  for (let i = 0; i < caseData.length; i++) {
    row = caseData[i];
    groupKeys = getGroupKeys(row, groupKey, dnaOrAa, coordinateMode);

    // If groupKeys is empty, that means that it's a sequence
    // that had its SNPs filtered out. We'll replace it with an "empty"
    // placeholder object to group by
    if (groupKeys.length === 0) {
      groupKeys = [-1];
    }

    groupKeys.forEach((_groupKey) => {
      // Replace the integer SNP ID with the actual SNP string
      if (groupKey === 'snp' && dnaOrAa === 'dna') {
        _groupKey = intToDnaSnp(_groupKey).snp_str;
      } else if (groupKey === 'snp' && dnaOrAa === 'aa') {
        if (coordinateMode === 'gene') {
          _groupKey = intToGeneAaSnp(_groupKey).snp_str;
        } else if (coordinateMode === 'protein') {
          _groupKey = intToProteinAaSnp(_groupKey).snp_str;
        }
      }

      // Create an entry for the group, if it doesn't already exist
      if (!Object.prototype.hasOwnProperty.call(aggCaseData, _groupKey)) {
        aggCaseData[_groupKey] = {};
      }
      // Create a date key, if it doesn't already exist
      if (
        !Object.prototype.hasOwnProperty.call(
          aggCaseData[_groupKey],
          row.collection_date
        )
      ) {
        aggCaseData[_groupKey][row.collection_date] = 0;
      }
      // Add case count
      aggCaseData[_groupKey][row.collection_date] += 1;
    });
  }

  // console.log(aggCaseData);

  // Compute total counts per group
  const countsPerGroup = {};
  Object.keys(aggCaseData).forEach((group) => {
    countsPerGroup[group] = Object.values(aggCaseData[group]).reduce(
      (memo, val) => memo + val,
      0
    );
  });

  // Expand obj of objs back into a list of objects
  const aggCaseDataList = [];
  Object.keys(aggCaseData).forEach((group) => {
    const dates = Object.keys(aggCaseData[group]);
    dates.forEach((date) => {
      aggCaseDataList.push({
        group: group,
        date: parseInt(date),
        cases_sum: aggCaseData[group][date],
        group_counts: countsPerGroup[group],
        color: getColorMethod(group),
      });
    });
  });

  // console.log(aggCaseDataList);

  return {
    aggCaseDataList,
    numSequencesBeforeMetadataFiltering,
    selectedRows,
    metadataCounts,
    aggLocationDataList,
  };
}

// Collapse case data by the grouping key
// i.e., collapse the date field for cases, so we can display group-wise stats
// in the data table
function aggCaseDataByGroup({
  caseData,
  coordinateMode,
  coordinateRanges,
  selectedGene,
  selectedProtein,
  groupKey,
  dnaOrAa,
  dateRange,
}) {
  const lineageCountObj = {};
  caseData.forEach((row) => {
    if (lineageCountObj[row.group]) lineageCountObj[row.group] += row.cases_sum;
    else {
      lineageCountObj[row.group] = row.cases_sum;
    }
  });

  const MAX_LINEAGE_SIZE = 10;
  let groupsToKeepObj;
  console.log(
    'lineage count obj: ',
    lineageCountObj,
    Object.keys(lineageCountObj).length > MAX_LINEAGE_SIZE
  );
  if (Object.keys(lineageCountObj).length > MAX_LINEAGE_SIZE) {
    let lineageCountArr = Object.entries(lineageCountObj);

    // this will sort it so that 0 is the biggest
    lineageCountArr.sort((a, b) => {
      if (a[1] < b[1]) {
        return 1;
      }
      if (a[1] > b[1]) {
        return -1;
      } else {
        return 0;
      }
    });

    lineageCountArr = lineageCountArr.slice(0, MAX_LINEAGE_SIZE);
    console.log('lineage count arr', lineageCountArr);
    groupsToKeepObj = Object.fromEntries(lineageCountArr);
    console.log('groups to keep', groupsToKeepObj);
  }

  // Aggregate case data by clade only (no dates)
  let caseDataAggGroup = {};
  let totalCaseCount = 0;

  // Filter by date
  caseData = filterByDate(caseData, dateRange);

  caseData.forEach((row) => {
    if (!Object.prototype.hasOwnProperty.call(caseDataAggGroup, row.group)) {
      caseDataAggGroup[row.group] = {};
      caseDataAggGroup[row.group]['cases_sum'] = 0;
    }
    caseDataAggGroup[row.group]['cases_sum'] += row.cases_sum;
    // Add to total case count
    totalCaseCount += row.cases_sum;
  });

  // Calculate percentages for cases
  // TODO: calculate growth rates
  Object.keys(caseDataAggGroup).forEach((row) => {
    caseDataAggGroup[row]['cases_percent'] =
      caseDataAggGroup[row]['cases_sum'] / totalCaseCount;
  });

  // We need a list of 0-indexed positions for the data table
  let changingPositions = {};
  // If we grouped by lineages/clades, then we need to find what SNPs
  // the lineages/clades correspond to, and add their SNP positions
  let groupSnpFunc = null;
  if (groupKey === 'lineage' || groupKey === 'clade') {
    if (dnaOrAa === 'dna') {
      groupSnpFunc = getDnaSnpsFromGroup.bind(this, groupKey);
    } else {
      if (coordinateMode === 'gene') {
        groupSnpFunc = getGeneAaSnpsFromGroup.bind(this, groupKey);
      } else if (coordinateMode === 'protein') {
        groupSnpFunc = getProteinAaSnpsFromGroup.bind(this, groupKey);
      }
    }

    // For each lineage/clade, add its changing positions
    let inRange = false;
    Object.keys(caseDataAggGroup).forEach((group) => {
      let groupSnps = groupSnpFunc(group);
      if (dnaOrAa === 'dna') {
        groupSnps.forEach((snp) => {
          inRange = _.some(coordinateRanges, (range) => {
            return snp.pos >= range[0] && snp.pos <= range[1];
          });

          if (inRange) {
            // positions are 1-indexed in the input file
            if (
              !Object.prototype.hasOwnProperty.call(
                changingPositions,
                snp.pos - 1
              )
            ) {
              changingPositions[snp.pos - 1] = {};
              changingPositions[snp.pos - 1]['ref'] = snp.ref;
            }
          }
        });
      } else {
        // AA-mode
        groupSnps.forEach((snp) => {
          if (coordinateMode === 'gene') {
            inRange = snp.gene === selectedGene.gene;
          } else if (coordinateMode === 'protein') {
            inRange = snp.protein === selectedProtein.protein;
          }

          if (inRange) {
            // positions are 0-indexed in the input file
            if (
              !Object.prototype.hasOwnProperty.call(changingPositions, snp.pos)
            ) {
              changingPositions[snp.pos] = {};
              if (coordinateMode === 'gene') {
                changingPositions[snp.pos]['gene'] = snp.gene;
              } else if (coordinateMode === 'protein') {
                changingPositions[snp.pos]['protein'] = snp.protein;
              }
              changingPositions[snp.pos]['ref'] = snp.ref;
            }
          }
        });
      }
    });
  } else if (groupKey === 'snp') {
    // If we're grouping by SNP, then just take the position embedded in each SNP string
    let snp_str_split = [];
    let pos = -1;
    if (dnaOrAa === 'dna') {
      // For DNA SNPs, the position is the first chunk (pos|ref|alt)
      // The DNA SNPs are 1-indexed, so -1 to make it 0-indexed
      Object.keys(caseDataAggGroup).forEach((snp_str) => {
        // Skip if the snp_str is Reference
        if (snp_str === 'Reference') {
          return;
        }
        snp_str_split = snp_str.split('|');
        pos = parseInt(snp_str_split[0]) - 1;

        if (!Object.prototype.hasOwnProperty.call(changingPositions, pos)) {
          changingPositions[pos] = {};
          changingPositions[pos]['ref'] = snp_str_split[1];
        }
      });
    } else {
      // For AA SNPs, the position is the second chunk (gene|pos|ref|alt)
      // The AA SNPs are 1-indexed, so -1 to make it 0-indexed
      Object.keys(caseDataAggGroup).forEach((snp_str) => {
        // Skip if the snp_str is Reference
        if (snp_str === 'Reference') {
          return;
        }

        snp_str_split = snp_str.split('|');
        pos = parseInt(snp_str_split[1]) - 1;

        if (!Object.prototype.hasOwnProperty.call(changingPositions, pos)) {
          changingPositions[pos] = {};
          if (coordinateMode === 'gene') {
            changingPositions[pos]['gene'] = snp_str_split[0];
          } else if (coordinateMode === 'protein') {
            changingPositions[pos]['protein'] = snp_str_split[0];
          }
          changingPositions[pos]['ref'] = snp_str_split[2];
        }
      });
    }
  }

  // Sort keys in changingPositions
  let changingPositionsOrdered = {};
  Object.keys(changingPositions)
    .sort((a, b) => parseInt(a) - parseInt(b))
    .forEach(function (key) {
      changingPositionsOrdered[key] = changingPositions[key];
    });
  changingPositions = changingPositionsOrdered;

  console.log(Object.keys(changingPositions).length, 'changing positions');
  //console.log(caseDataAggGroup);

  // Add each changing position as a new field for each row in caseDataAggGroup
  let ref_base = '';
  let alt_base = '';
  let snpRow = null;

  // Add the reference sequence, if it hasn't been added yet
  if (!Object.prototype.hasOwnProperty.call(caseDataAggGroup, 'Reference')) {
    caseDataAggGroup['Reference'] = {
      cases_sum: NaN,
      cases_percent: NaN,
    };
  }

  // For SNP data, break out the (gene), position, ref, and alt
  if (groupKey === 'snp') {
    let row_split;
    Object.keys(caseDataAggGroup).forEach((row) => {
      // Skip reference
      if (row === 'Reference') {
        // Since we're going to display the gene or pos later
        // in the data table, put these here so the reference
        // row renders correctly
        if (dnaOrAa == 'dna') {
          caseDataAggGroup[row]['pos'] = 'Reference';
        } else {
          if (coordinateMode === 'gene') {
            caseDataAggGroup[row]['gene'] = 'Reference';
          } else if (coordinateMode === 'protein') {
            caseDataAggGroup[row]['protein'] = 'Reference';
          }
        }
        return;
      }

      row_split = row.split('|');
      if (dnaOrAa === 'dna') {
        // DNA SNP string: pos|ref|alt
        caseDataAggGroup[row]['pos'] = parseInt(row_split[0]);
        caseDataAggGroup[row]['ref'] = row_split[1];
        caseDataAggGroup[row]['alt'] = row_split[2];
      } else {
        // AA SNP string: gene|pos|ref|alt or protein|pos|ref|alt
        if (coordinateMode === 'gene') {
          caseDataAggGroup[row]['gene'] = row_split[0];
        } else if (coordinateMode === 'protein') {
          caseDataAggGroup[row]['protein'] = row_split[0];
        }
        caseDataAggGroup[row]['pos'] = parseInt(row_split[1]);
        caseDataAggGroup[row]['ref'] = row_split[2];
        caseDataAggGroup[row]['alt'] = row_split[3];
      }
    });
  }

  Object.keys(caseDataAggGroup).forEach((group) => {
    Object.keys(changingPositions).forEach((pos) => {
      ref_base = changingPositions[pos]['ref'];
      alt_base = ref_base; // Default to same as reference
      pos = parseInt(pos);

      // Ignore finding alternate bases for the ref sequence
      if (group === 'Reference') {
        alt_base = ref_base;
      }
      // If we grouped by lineage, use the lineage name
      // to find a potential SNP at this location
      else if (groupKey === 'lineage' || groupKey === 'clade') {
        // lineageSnpFunc is defined above
        let groupSnps = groupSnpFunc(group);

        if (dnaOrAa === 'dna') {
          // The DNA SNPs are 0-indexed, so +1 to make it 1-indexed
          snpRow = _.findWhere(groupSnps, { pos: pos + 1 });
          if (snpRow !== undefined) {
            alt_base = snpRow.alt;
          }
        } else {
          // The AA SNPs are already 0-indexed
          snpRow = _.findWhere(groupSnps, { pos: pos });
          if (snpRow !== undefined) {
            alt_base = snpRow.alt;
          }
        }
      }
      // If we grouped by SNP, then the ref and alt base information
      // is embedded into the SNP string
      else if (groupKey === 'snp') {
        // DNA SNP string: pos|ref|alt
        if (dnaOrAa === 'dna') {
          if (pos === caseDataAggGroup[group]['pos'] - 1) {
            alt_base = caseDataAggGroup[group]['alt'];
          }
        }
        // AA SNP string: gene|pos|ref|alt
        else {
          if (pos === caseDataAggGroup[group]['pos'] - 1) {
            alt_base = caseDataAggGroup[group]['alt'];
          }
        }
      }

      caseDataAggGroup[group]['pos_' + pos.toString()] = alt_base;
    });
  });

  let getColorMethod;
  if (groupKey === 'lineage') {
    getColorMethod = getLineageColor;
  } else if (groupKey === 'clade') {
    getColorMethod = getCladeColor;
  } else if (groupKey === 'snp') {
    getColorMethod = getSnpColor;
  }

  // Object -> List of records
  Object.keys(caseDataAggGroup).forEach((group) => {
    caseDataAggGroup[group]['group'] = group;
    caseDataAggGroup[group]['color'] = getColorMethod(group);
    const parentkey = getParent(group);
    if (caseDataAggGroup[parentkey] && parentkey !== group) {
      caseDataAggGroup[group].parent = parentkey;
    } else if (group !== 'Reference') {
      caseDataAggGroup[group].parent = 'Reference';
    }
    caseDataAggGroup[group].name = group;
    caseDataAggGroup[group].id = group;
  });
  caseDataAggGroup = Object.values(caseDataAggGroup);

  return {
    caseDataAggGroup: caseDataAggGroup,
    changingPositions: changingPositions,
    groupsToKeepObj,
  };
}

/*
  processEntropyData() {
    let clade_id = -1;
    let clade_positions = [];
    let entropy_data = [];
    // For each row in the case data
    for(let i = 0; i < case_data.length; i += 1) {
    //for(let i = 0; i < 10000; i += 1) {
      // Lookup the positions for this clade_id
      clade_id = case_data[i]['clade_id'];
      clade_positions = _.findWhere(clade_data, { index: clade_id })['pos'];
      for(let j = 0; j < clade_positions.length; j += 1) {
        entropy_data.push({
          location_id: case_data[i]['location_id'],
          clade_id: clade_id,
          date: case_data[i]['date'],
          position: clade_positions[j],
          cases: case_data[i]['cases']
        });
      }
    }

    return entropy_data
  }
  */

self.addEventListener(
  'message',
  function (e) {
    const data = JSON.parse(e.data);

    let result;
    if (data.type === 'aggCaseDataByGroup') {
      //console.log('into casedata for agg', data);
      result = aggCaseDataByGroup(data);
    } else if (data.type === 'processCaseData') {
      //console.log('into casedata for process', data);
      result = processCaseData(data);
    }
    self.postMessage(JSON.stringify(result));
  },
  false
);
