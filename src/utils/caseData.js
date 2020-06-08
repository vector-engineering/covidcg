import initialCaseData from '../../data/case_data.json';
import {
  int_to_dna_snp,
  int_to_aa_snp,
  dna_snp_in_gene,
  aa_snp_in_gene,
} from './snpData';
import { loadLineageDnaSnp, loadLineageAaSnp } from './lineageData';
import _ from 'underscore';

const processedCaseData = _.map(initialCaseData, (row) => {
  row.date = new Date(row.date).getTime();
  return row;
});
// const initialLineageData = loadLineageData();

export function loadCaseData() {
  return processedCaseData;
}

function filterByLocation(caseData, locationIds) {
  // Skip if no locations selected
  if (locationIds.length === 0) {
    return caseData;
  }

  return _.filter(caseData, (row) => {
    return locationIds.includes(row.location_id);
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

  console.log('selected gene:', selectedGene);

  let newCaseData = [];

  if (groupKey === 'snp') {
    if (dnaOrAa === 'dna') {
      caseData.forEach((row) => {
        // Only keep SNPs that are within
        row['dna_snp_str'] = _.filter(row['dna_snp_str'], (snp_id) => {
          return dna_snp_in_gene(snp_id, selectedGene.gene);
        });
        newCaseData.push(row);
      });
    } else {
      caseData.forEach((row) => {
        // Only keep SNPs that are within
        row['aa_snp_str'] = _.filter(row['aa_snp_str'], (snp_id) => {
          return aa_snp_in_gene(snp_id, selectedGene.gene);
        });
        newCaseData.push(row);
      });
    }
  } /*else if (groupKey == 'snp_sig') {
  }*/

  return newCaseData;
}

export function processCaseData(locationIds, selectedGene, groupKey, dnaOrAa) {
  // let caseData = _.map(_caseData, (row) => Object.assign({}, row));
  let caseData = JSON.parse(JSON.stringify(processedCaseData));

  console.log('filtering by locationIds', locationIds);

  // Filter by location
  caseData = filterByLocation(caseData, locationIds);
  console.log(caseData.length, 'rows remaining after location filtering');
  // Filter by gene
  caseData = filterByGene(caseData, selectedGene, groupKey, dnaOrAa);
  console.log(caseData.length, 'rows remaining after gene filtering');

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
        _groupKey = int_to_dna_snp(_groupKey).snp_str;
      } else if (groupKey === 'snp' && dnaOrAa === 'aa') {
        _groupKey = int_to_aa_snp(_groupKey).snp_str;
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

  //console.log(aggCaseData);

  // Expand obj of objs back into a list of objects
  let aggCaseDataList = [];
  Object.keys(aggCaseData).forEach((group) => {
    let dates = Object.keys(aggCaseData[group]);
    dates.forEach((date) => {
      aggCaseDataList.push({
        group: group,
        date: date,
        cases_sum: aggCaseData[group][date],
      });
    });
  });

  //console.log(aggCaseDataList);

  return aggCaseDataList;
}

// Collapse case data by the grouping key
// i.e., collapse the date field for cases, so we can display group-wise stats
// in the data table
export function aggCaseDataByGroup(
  caseData,
  selectedGene,
  groupKey,
  dnaOrAa,
  dateRange
) {
  // Aggregate case data by clade only (no dates)
  let caseDataAggGroup = {};
  let totalCaseCount = 0;

  // Filter by date
  if (dateRange[0] > -1 && dateRange[1] > -1) {
    caseData = _.filter(caseData, (row) => {
      return row.date >= dateRange[0] && row.date <= dateRange[1];
    });
  }

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

  //console.log(caseData);
  //console.log(caseDataAggGroup);

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
        (row) => row.lineage == lineage
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
        // Skip if the snp_str is None
        if (snp_str === 'None') {
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
        // Skip if the snp_str is None
        if (snp_str === 'None') {
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

  console.log('Changing positions:', changingPositions);
  //console.log(caseDataAggGroup);

  // Add each changing position as a new field for each row in caseDataAggGroup
  let ref_base = '';
  let alt_base = '';
  let snpRow = null;

  // Add the reference sequence, if it hasn't been added yet
  if (!Object.prototype.hasOwnProperty.call(caseDataAggGroup, 'reference')) {
    caseDataAggGroup['Reference'] = {
      cases_sum: NaN,
      cases_percent: NaN,
    };
  }

  Object.keys(caseDataAggGroup).forEach((row) => {
    Object.keys(changingPositions).forEach((pos) => {
      ref_base = changingPositions[pos]['ref'];
      alt_base = ref_base; // Default to same as reference
      pos = parseInt(pos);

      // Ignore finding alternate bases for the ref sequence
      if (row === 'reference') {
        alt_base = ref_base;
      }
      // If we grouped by lineage, use the lineage name
      // to find a potential SNP at this location
      else if (groupKey === 'lineage') {
        lineageSnpData =
          dnaOrAa === 'dna' ? loadLineageDnaSnp() : loadLineageAaSnp();
        if (dnaOrAa === 'dna') {
          // The DNA SNPs are 1-indexed, so +1 to make it 1-indexed
          snpRow = _.findWhere(lineageSnpData, { lineage: row, pos: pos + 1 });
          if (snpRow !== undefined) {
            alt_base = snpRow.alt;
          }
        } else {
          // The AA SNPs are 0-indexed
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
          if (pos === parseInt(row.split('|')[0]) - 1) {
            alt_base = row.split('|')[2];
          }
        }
        // AA SNP string: gene|pos|ref|alt
        else {
          if (pos === parseInt(row.split('|')[1]) - 1) {
            alt_base = row.split('|')[3];
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

  // console.log(caseDataAggGroup);

  return {
    caseDataAggGroup: caseDataAggGroup,
    changingPositions: changingPositions,
  };

  /*

  // Obj to list of rows
  let caseDataAggGroupList = [];
  Object.keys(caseDataAggGroup).forEach((lineage) => {
    let lineageObj = {
      lineage: lineage === 'root' ? 'Reference' : lineage,
      cases_sum: caseDataAggGroup[lineage],
      cases_percent: caseDataAggGroup[lineage] / totalCaseCount,
      jmol: Math.random(),
    };

    let lineage_dat = _.filter(
      initialLineageData,
      (row) => row.lineage == lineage
    );

    changingPositions.forEach((pos) => {
      let ref_base = reference_seq[pos];
      let alt_base = ref_base;

      // Find the position in the SNPs for this lineage
      let row = _.findWhere(lineage_dat, { pos: pos + 1 });
      if (row !== undefined) {
        alt_base = row.alt;
      }

      lineageObj['pos_' + pos.toString()] = alt_base;
    });

    caseDataAggGroupList.push(lineageObj);
  });
  //console.log(caseDataAggGroupList);

  return {
    caseDataAggGroupList: caseDataAggGroupList,
    changingPositions: changingPositions,
  };*/
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
