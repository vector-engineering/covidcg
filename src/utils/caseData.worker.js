import initialCaseData from '../../data/case_data.json';
import { intToDnaSnp, intToAaSnp, dnaSnpInGene, aaSnpInGene } from './snpData';
import { loadLineageDnaSnp, loadLineageAaSnp } from './lineageData';
import {
  warmColors,
  coolColors,
  snpColorArray,
  incrementColor,
} from '../constants/colors';
import _ from 'underscore';

const processedCaseData = _.map(initialCaseData, (row) => {
  row.sample_date = new Date(row.sample_date).getTime();

  // Floor at 2020-01-01
  if (row.sample_date < 1577836800000) {
    row.sample_date = 1577836800000;
  }

  return row;
});

// function loadCaseData() {
//   return processedCaseData;
// }

let _warmColors = Object.assign({}, warmColors);
let _coolColors = Object.assign({}, coolColors);
let _snpColorArray = [...snpColorArray];

function iterateAndGetColor(colorobj, iter, i = 0) {
  if (i < iter) {
    return iterateAndGetColor(colorobj.child, iter, i + 1);
  } else {
    colorobj.colors.push(incrementColor(colorobj.colors.shift(), 1));
    return colorobj.colors[0];
  }
}

const getColor = _.memoize((group) => {
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

  const locationObj = convertToObj(locationIds);

  return caseData.filter((row) => {
    return locationObj[row.location_id] === 1;
  });
}

function filterByGene(caseData, selectedGene, groupKey, dnaOrAa) {
  // Don't need to do this if we're grouping by lineage
  if (groupKey === 'lineage') {
    return caseData;
  }

  // Don't need to do this if we selected all genes
  if (selectedGene.gene === 'All Genes') {
    return caseData;
  }

  console.log('selected gene:', selectedGene.gene);

  let newCaseData = [];

  if (groupKey === 'snp') {
    if (dnaOrAa === 'dna') {
      caseData.forEach((row) => {
        // Only keep SNPs that are within
        row['dna_snp_str'] = _.filter(row['dna_snp_str'], (snp_id) => {
          return dnaSnpInGene(snp_id, selectedGene.gene);
        });
        newCaseData.push(row);
      });
    } else {
      caseData.forEach((row) => {
        // Only keep SNPs that are within
        row['aa_snp_str'] = _.filter(row['aa_snp_str'], (snp_id) => {
          return aaSnpInGene(snp_id, selectedGene.gene);
        });
        newCaseData.push(row);
      });
    }
  } /*else if (groupKey == 'snp_sig') {
  }*/

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

function processCaseData(locationIds, selectedGene, groupKey, dnaOrAa) {
  // let caseData = _.map(_caseData, (row) => Object.assign({}, row));
  let caseData = JSON.parse(JSON.stringify(processedCaseData));
  //console.log('filtering by locationIds', locationIds);

  // Filter by location
  caseData = filterByLocation(caseData, locationIds);
  console.log(caseData.length, 'rows remaining after location filtering');
  // Filter by gene
  caseData = filterByGene(caseData, selectedGene, groupKey, dnaOrAa);
  console.log(caseData.length, 'rows remaining after gene filtering');

  // Get a list of Accession IDs and sample dates that are currently selected
  let selectedRows = _.map(caseData, (row) => {
    return row;
  });

  // Group by grouping key and sample date
  console.log('Grouping by', groupKey, 'and sample date');
  let aggCaseData = {};
  let row = {};
  let groupKeys = [];
  for (let i = 0; i < caseData.length; i++) {
    row = caseData[i];

    // If we're grouping by lineage or SNP signature, then the group
    // keys are literal
    if (groupKey === 'lineage') {
      groupKeys = [row['lineage']];
    } else if (groupKey === 'snp_sig') {
      if (dnaOrAa === 'dna') {
        groupKeys = [row['dna_snp_sig']];
      } else {
        groupKeys = [row['aa_snp_sig']];
      }
    }
    // If we're grouping by SNP, then the value is a semicolon-delimited
    // list of SNPs that we should treat separately
    else {
      if (dnaOrAa === 'dna') {
        groupKeys = row['dna_snp_str'];
      } else {
        groupKeys = row['aa_snp_str'];
      }
    }

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
        _groupKey = intToAaSnp(_groupKey).snp_str;
      }

      // Create an entry for the group, if it doesn't already exist
      if (!Object.prototype.hasOwnProperty.call(aggCaseData, _groupKey)) {
        aggCaseData[_groupKey] = {};
      }
      // Create a date key, if it doesn't already exist
      if (
        !Object.prototype.hasOwnProperty.call(
          aggCaseData[_groupKey],
          row.sample_date
        )
      ) {
        aggCaseData[_groupKey][row.sample_date] = 0;
      }
      // Add case count
      aggCaseData[_groupKey][row.sample_date] += 1;
    });
  }

  let getColorMethod = getColor;

  if (groupKey === 'snp') {
    getColorMethod = getSnpColor;
  }

  // Expand obj of objs back into a list of objects
  let aggCaseDataList = [];
  Object.keys(aggCaseData).forEach((group) => {
    let dates = Object.keys(aggCaseData[group]);
    dates.forEach((date) => {
      const color = getColorMethod(group);
      aggCaseDataList.push({
        group: group,
        date: date,
        cases_sum: aggCaseData[group][date],
        color,
      });
    });
  });

  //console.log(aggCaseDataList);

  return {
    aggCaseDataList: aggCaseDataList,
    selectedRows: selectedRows,
  };
}

// Collapse case data by the grouping key
// i.e., collapse the date field for cases, so we can display group-wise stats
// in the data table
function aggCaseDataByGroup({
  caseData,
  selectedGene,
  groupKey,
  dnaOrAa,
  dateRange,
}) {
  // Aggregate case data by clade only (no dates)
  let caseDataAggGroup = {};
  let totalCaseCount = 0;

  // Filter by date
  caseData = filterByDate(caseData, dateRange);
  console.log(caseData.length, 'rows remaining after date filtering');

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
  // If we grouped by lineages, then we need to find what SNPs
  // the lineages correspond to, and add their SNP positions
  let lineageSnpData;
  if (groupKey === 'lineage') {
    lineageSnpData =
      dnaOrAa === 'dna' ? loadLineageDnaSnp() : loadLineageAaSnp();

    // For each lineage, add its changing positions
    Object.keys(caseDataAggGroup).forEach((lineage) => {
      let lineageSnps = _.filter(
        lineageSnpData,
        (row) => row.lineage == lineage && row.group != 'Reference'
      );

      if (dnaOrAa === 'dna') {
        lineageSnps.forEach((row) => {
          if (row.pos > selectedGene.start && row.pos < selectedGene.end) {
            // positions are 1-indexed in the input file
            if (
              !Object.prototype.hasOwnProperty.call(
                changingPositions,
                row.pos - 1
              )
            ) {
              changingPositions[row.pos - 1] = {};
              changingPositions[row.pos - 1]['ref'] = row.ref;
            }
          }
        });
      } else {
        lineageSnps.forEach((row) => {
          if (row.gene === selectedGene.gene) {
            // positions are 0-indexed in the input file
            if (
              !Object.prototype.hasOwnProperty.call(changingPositions, row.pos)
            ) {
              changingPositions[row.pos] = {};
              changingPositions[row.pos]['gene'] = row.gene;
              changingPositions[row.pos]['ref'] = row.ref;
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
          changingPositions[pos]['gene'] = snp_str_split[0];
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
          caseDataAggGroup[row]['gene'] = 'Reference';
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
        // AA SNP string: gene|pos|ref|alt
        caseDataAggGroup[row]['gene'] = row_split[0];
        caseDataAggGroup[row]['pos'] = parseInt(row_split[1]);
        caseDataAggGroup[row]['ref'] = row_split[2];
        caseDataAggGroup[row]['alt'] = row_split[3];
      }
    });
  }

  Object.keys(caseDataAggGroup).forEach((row) => {
    Object.keys(changingPositions).forEach((pos) => {
      ref_base = changingPositions[pos]['ref'];
      alt_base = ref_base; // Default to same as reference
      pos = parseInt(pos);

      // Ignore finding alternate bases for the ref sequence
      if (row === 'Reference') {
        alt_base = ref_base;
      }
      // If we grouped by lineage, use the lineage name
      // to find a potential SNP at this location
      else if (groupKey === 'lineage') {
        lineageSnpData =
          dnaOrAa === 'dna' ? loadLineageDnaSnp() : loadLineageAaSnp();
        if (dnaOrAa === 'dna') {
          // The DNA SNPs are 0-indexed, so +1 to make it 1-indexed
          snpRow = _.findWhere(lineageSnpData, { lineage: row, pos: pos + 1 });
          if (snpRow !== undefined) {
            alt_base = snpRow.alt;
          }
        } else {
          // The AA SNPs are already 0-indexed
          snpRow = _.findWhere(lineageSnpData, { lineage: row, pos: pos });
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
          if (pos === caseDataAggGroup[row]['pos'] - 1) {
            alt_base = caseDataAggGroup[row]['alt'];
          }
        }
        // AA SNP string: gene|pos|ref|alt
        else {
          if (pos === caseDataAggGroup[row]['pos'] - 1) {
            alt_base = caseDataAggGroup[row]['alt'];
          }
        }
      }

      caseDataAggGroup[row]['pos_' + pos.toString()] = alt_base;
    });
  });

  // Object -> List of records
  Object.keys(caseDataAggGroup).forEach((group) => {
    caseDataAggGroup[group]['group'] = group;
  });
  caseDataAggGroup = Object.values(caseDataAggGroup);

  return {
    caseDataAggGroup: caseDataAggGroup,
    changingPositions: changingPositions,
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
      result = processCaseData(
        data.selectedLocationIds,
        data.selectedGene,
        data.groupKey,
        data.dnaOrAa
      );
    }
    self.postMessage(JSON.stringify(result));
  },
  false
);
