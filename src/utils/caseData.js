import initialCaseData from '../../processed_data/simulated_case_data.json';
import _ from 'underscore';

export function loadCaseData() {
  return initialCaseData;
}

export function processCaseData(caseData, cladeData, locationIds) {
  // Filter by location
  let filteredCaseData = caseData;
  if(locationIds.length > 0) {
    filteredCaseData = _.filter(filteredCaseData, row => {
      return locationIds.includes(row.location_id)
    });
  }

  // Group by clade_id and date, and sum up all cases
  let aggCaseData = {};
  let row = {};
  for(let i = 0; i < filteredCaseData.length; i++) {
    row = filteredCaseData[i];
    // Create a clade group, if it doesn't already exist
    if(!Object.prototype.hasOwnProperty.call(aggCaseData, row.clade_id)) {
      aggCaseData[row.clade_id] = {}
    }

    // Create a date key, if it doesn't already exist
    if(!Object.prototype.hasOwnProperty.call(aggCaseData[row.clade_id], row.date)) {
      aggCaseData[row.clade_id][row.date] = 0;
    } 
    // Add case count
    aggCaseData[row.clade_id][row.date] += row.cases;
  }
  // Expand obj of objs back into a list of objects
  let aggCaseDataList = [];
  Object.keys(aggCaseData).forEach(cladeId => {
    let dates = Object.keys(aggCaseData[cladeId]);
    dates.forEach(date => {
      aggCaseDataList.push({
        clade_name: _.findWhere(cladeData, { index: parseInt(cladeId) }).clade,
        date: date,
        cases_sum: aggCaseData[cladeId][date]
      });
    });
  });

  //console.log(aggCaseDataList);

  return(aggCaseDataList);
}

// Collapse case data by clade
export function aggCaseDataByClade(caseData, cladeData, start_pos, end_pos, dateRange) {
  // Aggregate case data by clade only (no dates)
  let caseDataAggClade = {};
  let totalCaseCount = 0;

  // Filter by date
  if(dateRange[0] > -1 && dateRange[1] > -1) {
    caseData = _.filter(caseData, row => {
      return row.date >= dateRange[0] && row.date <= dateRange[1]
    });
  }

  caseData.forEach(row => {
    if(!Object.prototype.hasOwnProperty.call(caseDataAggClade, row.clade_name)) {
      caseDataAggClade[row.clade_name] = 0;
    }
    caseDataAggClade[row.clade_name] += row.cases_sum;
    // Add to total case count
    totalCaseCount += row.cases_sum;
  });

  // Add the reference sequence, if it hasn't been added yet
  if(!Object.prototype.hasOwnProperty.call(caseDataAggClade, 'root')) {
    caseDataAggClade['root'] = 0;
  }

  // Get all changing positions
  let changingPositions = new Set();
  Object.keys(caseDataAggClade).forEach(clade_name => {
    let clade = _.findWhere(cladeData, { clade: clade_name });
    clade.pos.forEach(pos => {
      if(pos > start_pos && pos < end_pos) {
        // positions are 1-indexed in the input file
        changingPositions.add(pos - 1);
      }
    });
  });
  // Convert from set to sorted array
  changingPositions = Array.from(changingPositions).sort((a, b) => a - b);

  // Obj to list of rows
  let caseDataAggCladeList = [];
  Object.keys(caseDataAggClade).forEach(clade_name => {

    let cladeObj = {
      clade_name: clade_name === 'root' ? 'Reference' : clade_name,
      cases_sum: caseDataAggClade[clade_name],
      cases_percent: caseDataAggClade[clade_name] / totalCaseCount
    };

    // Get the bases for this clade's sequence at each position 
    // in changingPositions
    let cladeSeq = _.findWhere(cladeData, { clade: clade_name })['seq'];
    // For each position, add a key to the table object
    changingPositions.forEach(pos => {
      cladeObj['pos_' + pos.toString()] = cladeSeq[pos];
    });

    caseDataAggCladeList.push(cladeObj);
  });
  //console.log(caseDataAggCladeList);

  return {
    caseDataAggCladeList: caseDataAggCladeList,
    changingPositions: changingPositions
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
