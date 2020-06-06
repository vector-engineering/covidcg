import initialCaseData from '../../data/case_data.json';
import refSeq from '../../data/reference.json';
import _ from 'underscore';

const reference_seq = refSeq['ref_seq'];

export function loadCaseData() {
  return initialCaseData;
}

export function processCaseData(caseData, locationIds) {
  // Filter by location
  let filteredCaseData = caseData;
  if (locationIds.length > 0) {
    filteredCaseData = _.filter(filteredCaseData, (row) => {
      return locationIds.includes(row.loc_id);
    });
  }

  // Group by lineage and date, and sum up all cases
  let aggCaseData = {};
  let row = {};
  for (let i = 0; i < filteredCaseData.length; i++) {
    row = filteredCaseData[i];
    // Create a lineage group, if it doesn't already exist
    if (!Object.prototype.hasOwnProperty.call(aggCaseData, row.lineage)) {
      aggCaseData[row.lineage] = {};
    }

    // Create a date key, if it doesn't already exist
    if (
      !Object.prototype.hasOwnProperty.call(aggCaseData[row.lineage], row.date)
    ) {
      aggCaseData[row.lineage][row.date] = 0;
    }
    // Add case count
    aggCaseData[row.lineage][row.date] += row.cases;
  }
  // Expand obj of objs back into a list of objects
  let aggCaseDataList = [];
  Object.keys(aggCaseData).forEach((lineage) => {
    let dates = Object.keys(aggCaseData[lineage]);
    dates.forEach((date) => {
      aggCaseDataList.push({
        lineage: lineage,
        date: date,
        cases_sum: aggCaseData[lineage][date],
      });
    });
  });

  //console.log(aggCaseDataList);

  return aggCaseDataList;
}

// Collapse case data by clade
export function aggCaseDataByLineage(
  caseData,
  lineageData,
  start_pos,
  end_pos,
  dateRange
) {
  // Aggregate case data by clade only (no dates)
  let caseDataAggLineage = {};
  let totalCaseCount = 0;

  // Filter by date
  if (dateRange[0] > -1 && dateRange[1] > -1) {
    caseData = _.filter(caseData, (row) => {
      return row.date >= dateRange[0] && row.date <= dateRange[1];
    });
  }

  caseData.forEach((row) => {
    if (
      !Object.prototype.hasOwnProperty.call(caseDataAggLineage, row.lineage)
    ) {
      caseDataAggLineage[row.lineage] = 0;
    }
    caseDataAggLineage[row.lineage] += row.cases_sum;
    // Add to total case count
    totalCaseCount += row.cases_sum;
  });

  // Add the reference sequence, if it hasn't been added yet
  if (!Object.prototype.hasOwnProperty.call(caseDataAggLineage, 'root')) {
    caseDataAggLineage['root'] = 0;
  }

  // Get all changing positions
  let changingPositions = new Set();
  Object.keys(caseDataAggLineage).forEach((lineage) => {
    let lineage_dat = _.filter(lineageData, (row) => row.lineage == lineage);

    lineage_dat.forEach((row) => {
      if (row.pos > start_pos && row.pos < end_pos) {
        // positions are 1-indexed in the input file
        changingPositions.add(row.pos - 1);
      }
    });
  });
  // Convert from set to sorted array
  changingPositions = Array.from(changingPositions).sort((a, b) => a - b);

  // Obj to list of rows
  let caseDataAggLineageList = [];
  Object.keys(caseDataAggLineage).forEach((lineage) => {
    let lineageObj = {
      lineage: lineage === 'root' ? 'Reference' : lineage,
      cases_sum: caseDataAggLineage[lineage],
      cases_percent: caseDataAggLineage[lineage] / totalCaseCount,
      jmol: Math.random(),
    };

    let lineage_dat = _.filter(lineageData, (row) => row.lineage == lineage);

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

    caseDataAggLineageList.push(lineageObj);
  });
  //console.log(caseDataAggLineageList);

  return {
    caseDataAggLineageList: caseDataAggLineageList,
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
