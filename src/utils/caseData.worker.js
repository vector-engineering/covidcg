import initialCaseData from '../../data/case_data.json';
import {
  intToDnaSnp,
  intToGeneAaSnp,
  intToProteinAaSnp,
  getSnvColor,
  formatSnv,
} from './snpData';
import {
  getDnaSnpsFromGroup,
  getGeneAaSnpsFromGroup,
  getProteinAaSnpsFromGroup,
  getLineageColor,
  getCladeColor,
} from './lineageData';
import { getLocationIds } from './location';
import { dataDate } from './version';
import { aggregate } from './transform';

import { countMetadataFields } from './metadata';
import _ from 'underscore';

import { GROUP_KEYS, DNA_OR_AA, COORDINATE_MODES } from '../constants/config';
import { REFERENCE_GROUP } from '../constants/groups';

const dataDateInt = new Date(dataDate).getTime();
const processedCaseData = _.reject(
  _.map(initialCaseData, (row) => {
    row.collection_date = new Date(row.collection_date).getTime();
    return row;
  }),
  (row) => {
    // Remove cases before 2019-12-15 and after the dataDate
    return (
      row.collection_date < 1576368000000 || row.collection_date > dataDateInt
    );
  }
);

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

  if (coordinateMode === COORDINATE_MODES.COORD_GENE) {
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
  } else if (coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
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
function filterByDate(caseData, dateRange, dateKey = 'date') {
  // Filter by date
  if (dateRange[0] > -1 && dateRange[1] > -1) {
    return _.filter(caseData, (row) => {
      return row[dateKey] >= dateRange[0] && row[dateKey] <= dateRange[1];
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
  if (groupKey === GROUP_KEYS.GROUP_LINEAGE) {
    return [row['lineage']];
  } else if (groupKey === GROUP_KEYS.GROUP_CLADE) {
    return [row['clade']];
  }
  // If we're grouping by SNP, then the value is a semicolon-delimited
  // list of SNPs that we should treat separately
  else if (groupKey === GROUP_KEYS.GROUP_SNV) {
    if (dnaOrAa === DNA_OR_AA.DNA) {
      return row['dna_snp_str'];
    } else {
      if (coordinateMode === COORDINATE_MODES.COORD_GENE) {
        return row['gene_aa_snp_str'];
      } else if (coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
        return row['protein_aa_snp_str'];
      }
    }
  }
}

function processCaseData({
  selectedLocationNodes,
  coordinateMode,
  coordinateRanges,
  selectedGene,
  selectedProtein,
  groupKey,
  dnaOrAa,
  selectedMetadataFields,
  ageRange,
  dateRange,
}) {
  // let caseData = _.map(_caseData, (row) => Object.assign({}, row));
  let caseData = JSON.parse(JSON.stringify(processedCaseData));

  // Filter by location
  const selectedLocationIds = getLocationIds(selectedLocationNodes);
  caseData = filterByLocation(caseData, selectedLocationIds);
  // console.log(caseData.length, 'rows remaining after location filtering');
  // Filter by coordinate range (DNA mode only)
  if (groupKey === GROUP_KEYS.GROUP_SNV && dnaOrAa === DNA_OR_AA.DNA) {
    caseData = filterByCoordinateRange({ caseData, coordinateRanges });
    // console.log(
    //   caseData.length,
    //   'rows remaining after coordinate range filtering'
    // );
  }
  // Filter by gene/protein (AA mode only, gene or protein selected)
  if (groupKey === GROUP_KEYS.GROUP_SNV && dnaOrAa === DNA_OR_AA.AA) {
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

  const metadataCountsAfterFiltering = countMetadataFields(caseData);

  // For the location tab:
  const aggCaseData = {};
  const countsPerLocation = {};
  const countsPerLocationDate = {};

  // Build a map of location_id --> node
  // Also while we're at it, create an entry for this node
  // in our data objects
  const locationIdToNodeMap = {};
  for (let i = 0; i < selectedLocationNodes.length; i++) {
    selectedLocationIds[i].forEach((locationId) => {
      locationIdToNodeMap[locationId] = selectedLocationNodes[i].value;
    });
    countsPerLocation[selectedLocationNodes[i].value] = 0;
    countsPerLocationDate[selectedLocationNodes[i].value] = {};
  }

  // If we have selected groups (lineages/snps/clades), then filter for that
  const uniqueGroupKeys = new Set();
  let groupKeys = [];
  let location;
  caseData.forEach((row) => {
    countsPerLocation[locationIdToNodeMap[row.location_id]] += 1;
    !(
      row.collection_date in
      countsPerLocationDate[locationIdToNodeMap[row.location_id]]
    ) &&
      (countsPerLocationDate[locationIdToNodeMap[row.location_id]][
        row.collection_date
      ] = 0);
    countsPerLocationDate[locationIdToNodeMap[row.location_id]][
      row.collection_date
    ] += 1;
    groupKeys = getGroupKeys(row, groupKey, dnaOrAa, coordinateMode);

    // If groupKeys is empty, that means that it's a sequence
    // that had its SNPs filtered out. We'll replace it with an "empty"
    // placeholder object to group by
    if (groupKeys.length === 0) {
      groupKeys = [-1];
    }

    location = locationIdToNodeMap[row.location_id];
    !(location in aggCaseData) && (aggCaseData[location] = {});
    !(row.collection_date in aggCaseData[location]) &&
      (aggCaseData[location][row.collection_date] = {});
    groupKeys.forEach((group) => {
      // Replace the integer SNP ID with the actual SNP string
      if (groupKey === GROUP_KEYS.GROUP_SNV && dnaOrAa === DNA_OR_AA.DNA) {
        group = intToDnaSnp(group).snp_str;
      } else if (
        groupKey === GROUP_KEYS.GROUP_SNV &&
        dnaOrAa === DNA_OR_AA.AA
      ) {
        if (coordinateMode === COORDINATE_MODES.COORD_GENE) {
          group = intToGeneAaSnp(group).snp_str;
        } else if (coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
          group = intToProteinAaSnp(group).snp_str;
        }
      }

      !(group in aggCaseData[location][row.collection_date]) &&
        (aggCaseData[location][row.collection_date][group] = 0);
      aggCaseData[location][row.collection_date][group] += 1;

      uniqueGroupKeys.add(group);
    });
  });

  // Process the object of counts per location, per date
  // const countsPerLocationDateList = [];
  // Object.keys(countsPerLocationDate).forEach((location) => {
  //   Object.keys(countsPerLocationDate[location]).forEach((date) => {
  //     countsPerLocationDateList.push({
  //       location: location,
  //       date: parseInt(date),
  //       count: countsPerLocationDate[location][date],
  //     });
  //   });
  // });

  let getColorMethod;
  if (groupKey === GROUP_KEYS.GROUP_LINEAGE) {
    getColorMethod = getLineageColor;
  } else if (groupKey === GROUP_KEYS.GROUP_CLADE) {
    getColorMethod = getCladeColor;
  } else if (groupKey === GROUP_KEYS.GROUP_SNV) {
    getColorMethod = getSnvColor;
  }

  const dataAggLocationGroupDate = [];
  Object.keys(aggCaseData).forEach((location) => {
    const dates = Object.keys(aggCaseData[location]);
    dates.forEach((date) => {
      groupKeys = Object.keys(aggCaseData[location][date]);
      groupKeys.forEach((group) => {
        dataAggLocationGroupDate.push({
          location: location,
          date: parseInt(date),
          group: group,
          groupName:
            groupKey === GROUP_KEYS.GROUP_SNV
              ? formatSnv(group, dnaOrAa)
              : group,
          cases_sum: aggCaseData[location][date][group],
          location_counts: countsPerLocation[location],
          color: getColorMethod(group),
        });
      });
    });
  });

  // Get a list of Accession IDs and sample dates that are currently selected
  const selectedAccessionIds = [];
  const selectedAckIds = [];
  caseData = filterByDate(caseData, dateRange, 'collection_date');
  caseData.forEach((row) => {
    selectedAccessionIds.push(row['Accession ID']);
    selectedAckIds.push(row['ack_id']);
  });

  // Aggregate by group and date
  const dataAggGroupDate = aggregate({
    data: dataAggLocationGroupDate,
    groupby: ['date', 'group', 'groupName'],
    fields: ['cases_sum', 'color', 'location_counts'],
    ops: ['sum', 'max', 'max'],
    as: ['cases_sum', 'color', 'location_counts'],
  });

  return {
    filteredCaseData: caseData,
    dataAggLocationGroupDate,
    dataAggGroupDate,
    numSequencesBeforeMetadataFiltering,
    metadataCounts,
    metadataCountsAfterFiltering,
    selectedAccessionIds,
    selectedAckIds,

    countsPerLocation,
    //countsPerLocationDate: countsPerLocationDateList,
    countsPerLocationDate,
  };
}

// Collapse case data by the grouping key
// i.e., collapse the date field for cases, so we can display group-wise stats
// in the data table
function aggCaseDataByGroup({
  totalSequenceCount,
  dataAggGroupDate,
  coordinateMode,
  coordinateRanges,
  selectedGene,
  selectedProtein,
  groupKey,
  dnaOrAa,
  dateRange,
}) {
  // console.log(dateRange);

  let getColorMethod;
  if (groupKey === GROUP_KEYS.GROUP_LINEAGE) {
    getColorMethod = getLineageColor;
  } else if (groupKey === GROUP_KEYS.GROUP_CLADE) {
    getColorMethod = getCladeColor;
  } else if (groupKey === GROUP_KEYS.GROUP_SNV) {
    getColorMethod = getSnvColor;
  }

  const groupCountObj = {};
  dataAggGroupDate.forEach((row) => {
    if (groupCountObj[row.group]) groupCountObj[row.group] += row.cases_sum;
    else groupCountObj[row.group] = row.cases_sum;
  });

  let groupCountArr = Object.entries(groupCountObj);
  groupCountArr.forEach((item) => {
    item.push(getColorMethod(item[0]));
    // Push the formatted group/SNV
    item.push(
      groupKey === GROUP_KEYS.GROUP_SNV ? formatSnv(item[0], dnaOrAa) : item[0]
    );
  });

  // this will sort it so that 0 is the biggest
  groupCountArr.sort((a, b) => {
    if (a[1] < b[1]) {
      return 1;
    }
    if (a[1] > b[1]) {
      return -1;
    } else {
      return 0;
    }
  });

  // Filter by date
  dataAggGroupDate = filterByDate(dataAggGroupDate, dateRange);

  // Create the same group count object, but after date filtering
  const groupCountDateFilteredObj = {};
  dataAggGroupDate.forEach((row) => {
    if (groupCountDateFilteredObj[row.group])
      groupCountDateFilteredObj[row.group] += row.cases_sum;
    else groupCountDateFilteredObj[row.group] = row.cases_sum;
  });
  let groupCountDateFilteredArr = Object.entries(groupCountDateFilteredObj);
  groupCountDateFilteredArr.forEach((item) => {
    item.push(getColorMethod(item[0]));
    // Push the formatted group/SNV
    item.push(
      groupKey === GROUP_KEYS.GROUP_SNV ? formatSnv(item[0], dnaOrAa) : item[0]
    );
  });
  groupCountDateFilteredArr.sort((a, b) => {
    if (a[1] < b[1]) {
      return 1;
    }
    if (a[1] > b[1]) {
      return -1;
    } else {
      return 0;
    }
  });

  // Aggregate case data by clade only (no dates)
  let dataAggGroup = {};
  dataAggGroupDate.forEach((row) => {
    if (!Object.prototype.hasOwnProperty.call(dataAggGroup, row.group)) {
      dataAggGroup[row.group] = {};
      dataAggGroup[row.group]['cases_sum'] = 0;
    }
    dataAggGroup[row.group]['cases_sum'] += row.cases_sum;
  });

  // Calculate percentages for cases
  // TODO: calculate growth rates
  Object.keys(dataAggGroup).forEach((row) => {
    dataAggGroup[row]['cases_percent'] =
      dataAggGroup[row]['cases_sum'] / totalSequenceCount;
  });

  // We need a list of 0-indexed positions for the data table
  let changingPositions = {};
  // If we grouped by lineages/clades, then we need to find what SNPs
  // the lineages/clades correspond to, and add their SNP positions
  let groupSnpFunc = null;
  if (
    groupKey === GROUP_KEYS.GROUP_LINEAGE ||
    groupKey === GROUP_KEYS.GROUP_CLADE
  ) {
    if (dnaOrAa === DNA_OR_AA.DNA) {
      groupSnpFunc = getDnaSnpsFromGroup.bind(this, groupKey);
    } else {
      if (coordinateMode === COORDINATE_MODES.COORD_GENE) {
        groupSnpFunc = getGeneAaSnpsFromGroup.bind(this, groupKey);
      } else if (coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
        groupSnpFunc = getProteinAaSnpsFromGroup.bind(this, groupKey);
      }
    }

    // For each lineage/clade, add its changing positions
    let inRange = false;
    Object.keys(dataAggGroup).forEach((group) => {
      let groupSnps = groupSnpFunc(group);
      if (dnaOrAa === DNA_OR_AA.DNA) {
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
          if (coordinateMode === COORDINATE_MODES.COORD_GENE) {
            inRange = snp.gene === selectedGene.gene;
          } else if (coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
            inRange = snp.protein === selectedProtein.protein;
          }

          if (inRange) {
            // positions are 0-indexed in the input file
            if (
              !Object.prototype.hasOwnProperty.call(changingPositions, snp.pos)
            ) {
              changingPositions[snp.pos] = {};
              if (coordinateMode === COORDINATE_MODES.COORD_GENE) {
                changingPositions[snp.pos]['gene'] = snp.gene;
              } else if (coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
                changingPositions[snp.pos]['protein'] = snp.protein;
              }
              changingPositions[snp.pos]['ref'] = snp.ref;
            }
          }
        });
      }
    });
  } else if (groupKey === GROUP_KEYS.GROUP_SNV) {
    // If we're grouping by SNP, then just take the position embedded in each SNP string
    let snp_str_split = [];
    let pos = -1;
    if (dnaOrAa === DNA_OR_AA.DNA) {
      // For DNA SNPs, the position is the first chunk (pos|ref|alt)
      // The DNA SNPs are 1-indexed, so -1 to make it 0-indexed
      Object.keys(dataAggGroup).forEach((snp_str) => {
        // Skip if the snp_str is Reference
        if (snp_str === REFERENCE_GROUP) {
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
      Object.keys(dataAggGroup).forEach((snp_str) => {
        // Skip if the snp_str is Reference
        if (snp_str === REFERENCE_GROUP) {
          return;
        }

        snp_str_split = snp_str.split('|');
        pos = parseInt(snp_str_split[1]) - 1;

        if (!Object.prototype.hasOwnProperty.call(changingPositions, pos)) {
          changingPositions[pos] = {};
          if (coordinateMode === COORDINATE_MODES.COORD_GENE) {
            changingPositions[pos]['gene'] = snp_str_split[0];
          } else if (coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
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

  // console.log(Object.keys(changingPositions).length, 'changing positions');
  //console.log(dataAggGroup);

  // Add each changing position as a new field for each row in dataAggGroup
  let ref_base = '';
  let alt_base = '';
  let snpRow = null;

  // Add the reference sequence, if it hasn't been added yet
  if (!Object.prototype.hasOwnProperty.call(dataAggGroup, REFERENCE_GROUP)) {
    dataAggGroup[REFERENCE_GROUP] = {
      cases_sum: NaN,
      cases_percent: NaN,
    };
  }

  // For SNP data, break out the (gene), position, ref, and alt
  if (groupKey === GROUP_KEYS.GROUP_SNV) {
    let row_split;
    Object.keys(dataAggGroup).forEach((row) => {
      // Skip reference
      if (row === REFERENCE_GROUP) {
        // Since we're going to display the gene or pos later
        // in the data table, put these here so the reference
        // row renders correctly
        if (dnaOrAa == DNA_OR_AA.DNA) {
          dataAggGroup[row]['pos'] = REFERENCE_GROUP;
        } else {
          if (coordinateMode === COORDINATE_MODES.COORD_GENE) {
            dataAggGroup[row]['gene'] = REFERENCE_GROUP;
          } else if (coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
            dataAggGroup[row]['protein'] = REFERENCE_GROUP;
          }
        }
        return;
      }

      row_split = row.split('|');
      if (dnaOrAa === DNA_OR_AA.DNA) {
        // DNA SNP string: pos|ref|alt
        dataAggGroup[row]['pos'] = parseInt(row_split[0]);
        dataAggGroup[row]['ref'] = row_split[1];
        dataAggGroup[row]['alt'] = row_split[2];
      } else {
        // AA SNP string: gene|pos|ref|alt or protein|pos|ref|alt
        if (coordinateMode === COORDINATE_MODES.COORD_GENE) {
          dataAggGroup[row]['gene'] = row_split[0];
        } else if (coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
          dataAggGroup[row]['protein'] = row_split[0];
        }
        dataAggGroup[row]['pos'] = parseInt(row_split[1]);
        dataAggGroup[row]['ref'] = row_split[2];
        dataAggGroup[row]['alt'] = row_split[3];
      }
    });
  }

  Object.keys(dataAggGroup).forEach((group) => {
    Object.keys(changingPositions).forEach((pos) => {
      ref_base = changingPositions[pos]['ref'];
      alt_base = ref_base; // Default to same as reference
      pos = parseInt(pos);

      // Ignore finding alternate bases for the ref sequence
      if (group === REFERENCE_GROUP) {
        alt_base = ref_base;
      }
      // If we grouped by lineage, use the lineage name
      // to find a potential SNP at this location
      else if (
        groupKey === GROUP_KEYS.GROUP_LINEAGE ||
        groupKey === GROUP_KEYS.GROUP_CLADE
      ) {
        // lineageSnpFunc is defined above
        let groupSnps = groupSnpFunc(group);

        if (dnaOrAa === DNA_OR_AA.DNA) {
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
      else if (groupKey === GROUP_KEYS.GROUP_SNV) {
        // DNA SNP string: pos|ref|alt
        if (dnaOrAa === DNA_OR_AA.DNA) {
          if (pos === dataAggGroup[group]['pos'] - 1) {
            alt_base = dataAggGroup[group]['alt'];
          }
        }
        // AA SNP string: gene|pos|ref|alt
        else {
          if (pos === dataAggGroup[group]['pos'] - 1) {
            alt_base = dataAggGroup[group]['alt'];
          }
        }
      }

      dataAggGroup[group]['pos_' + pos.toString()] = alt_base;
    });
  });

  // Object -> List of records
  Object.keys(dataAggGroup).forEach((group) => {
    dataAggGroup[group]['group'] = group;
    dataAggGroup[group]['groupName'] =
      groupKey === GROUP_KEYS.GROUP_SNV ? formatSnv(group, dnaOrAa) : group;
    dataAggGroup[group]['color'] = getColorMethod(group);
    const parentkey = getParent(group);
    if (dataAggGroup[parentkey] && parentkey !== group) {
      dataAggGroup[group].parent = parentkey;
    } else if (group !== REFERENCE_GROUP) {
      dataAggGroup[group].parent = REFERENCE_GROUP;
    }
    dataAggGroup[group].name = group;
    dataAggGroup[group].id = group;
  });
  dataAggGroup = Object.values(dataAggGroup);

  return {
    dataAggGroup: dataAggGroup,
    changingPositions: changingPositions,
    groupCountArr,
    groupCountDateFilteredArr,
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
    const data = e.data;

    let result;
    if (data.type === 'aggCaseDataByGroup') {
      //console.log('into casedata for agg', data);
      result = aggCaseDataByGroup(data);
    } else if (data.type === 'processCaseData') {
      //console.log('into casedata for process', data);
      result = processCaseData(data);
    }
    self.postMessage(result);
  },
  false
);
