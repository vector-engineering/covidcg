import {
  GROUP_SNV,
  DNA_OR_AA,
  COORDINATE_MODES,
  GROUPS,
} from '../../constants/defs.json';
import { config } from '../../config';

function downloadAggCaseDataGroup({
  groupKey,
  dataAggGroup,
  changingPositions,
  coordinateMode,

  // SNV data
  intToDnaSnvMap,
  intToGeneAaSnvMap,
  intToProteinAaSnvMap,

  // Lineage data
  groupSnvMap,
  // groupColorMap,
}) {
  let csvString = '';

  // Write headers
  csvString = 'lineage,seqs,seqs_percent,nt_snps,aa_snps,';
  // Add position column headers
  csvString += changingPositions.map((pos) => (pos + 1).toString()).join(',');
  csvString += '\n';

  for (let i = 0; i < dataAggGroup.length; i++) {
    let row = dataAggGroup[i];
    // Skip if it's the reference row
    if (row['group'] === GROUPS.REFERENCE_GROUP) {
      continue;
    }

    // Write lineage and counts
    csvString +=
      row['group'] + ',' + row['counts'] + ',' + row['percent'] + ',';

    // If it's the "Other" row, then don't try to get any SNVs
    if (row['group'] === GROUPS.OTHER_GROUP) {
      csvString += ',,';
      // Add empty field for each changing position
      csvString += changingPositions.map(() => '').join(',');
      csvString += '\n';
      continue;
    }

    // Get NT SNPs
    let ntSnps = groupSnvMap[groupKey][row['group']]['dna_snp_ids'].map(
      (snvId) => {
        return intToDnaSnvMap[snvId];
      }
    );
    // Skip if it's empty
    if (ntSnps.length === 0) {
      csvString += ',';
    } else {
      // Loop thru SNPs and print them as a list
      ntSnps = ntSnps.map((snp) => {
        // Format as pos|ref|alt
        // Position is 0-indexed, so make it 1-indexed
        return (
          (parseInt(snp['pos']) + 1).toString() +
          '|' +
          snp['ref'] +
          '|' +
          snp['alt']
        );
      });
      csvString += '"[' + ntSnps.join(',') + ']",';
    }

    // Get AA SNPs
    let aaSnps = [];
    if (coordinateMode === COORDINATE_MODES.COORD_GENE) {
      aaSnps = groupSnvMap[groupKey][row['group']]['gene_aa_snp_ids'].map(
        (snvId) => {
          return intToGeneAaSnvMap[snvId];
        }
      );
    } else if (coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
      aaSnps = groupSnvMap[groupKey][row['group']]['protein_aa_snp_ids'].map(
        (snvId) => {
          return intToProteinAaSnvMap[snvId];
        }
      );
    }
    // Skip if it's empty
    if (aaSnps.length === 0) {
      csvString += ',';
    } else {
      // Loop thru SNPs and print them as a list
      aaSnps = aaSnps.map((snp) => {
        // Format as gene|pos|ref|alt
        // Position is 0-indexed, so make it 1-indexed
        let label = '';
        if (coordinateMode === COORDINATE_MODES.COORD_GENE) {
          label = snp['gene'];
        } else if (coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
          label = snp['protein'];
        }
        return (
          label +
          '|' +
          parseInt(snp['pos']).toString() +
          '|' +
          snp['ref'] +
          '|' +
          snp['alt']
        );
      });
      csvString += '"[' + aaSnps.join(',') + ']",';
    }

    // Add letters at positions
    csvString += changingPositions
      .map((pos) => row['pos_' + pos.toString()])
      .join(',');

    csvString += '\n';
  }

  return csvString;
}

function downloadAggCaseDataSnp(dnaOrAa, dataAggGroup, changingPositions) {
  let csvString = '';

  // Write headers
  csvString = 'snp_str,';
  // Add DNA headers
  if (dnaOrAa === DNA_OR_AA.DNA) {
    csvString += 'pos,ref,alt,';
  }
  // Add AA headers - just DNA but with gene too
  else {
    csvString += 'gene,pos,ref,alt,';
  }
  // Add count headers
  csvString += 'seqs,seqs_percent,';
  // Add position column headers
  csvString += changingPositions.map((pos) => (pos + 1).toString()).join(',');
  csvString += '\n';

  for (let i = 0; i < dataAggGroup.length; i++) {
    let row = dataAggGroup[i];

    // Write the SNP string
    // For DNA, its pos|ref|alt
    // For AA, its gene|pos|ref|alt
    // And then write the SNP chunks
    if (dnaOrAa === DNA_OR_AA.DNA) {
      // Handle reference row
      if (Object.values(GROUPS).includes(row['group'])) {
        csvString += row['group'] + ',,,,';
      } else {
        csvString += [row['pos'], row['ref'], row['alt']].join('|') + ',';
        csvString += [row['pos'], row['ref'], row['alt']].join(',') + ',';
      }
    } else {
      // Handle reference row
      if (Object.values(GROUPS).includes(row['group'])) {
        csvString += row['group'] + ',,,,,';
      } else {
        csvString +=
          [row['gene'], row['pos'], row['ref'], row['alt']].join('|') + ',';
        csvString +=
          [row['gene'], row['pos'], row['ref'], row['alt']].join(',') + ',';
      }
    }

    // Write the sequence counts/percents
    csvString += row['counts'] + ',' + row['percent'] + ',';

    // Add letters at positions
    csvString += changingPositions
      .map((pos) => row['pos_' + pos.toString()])
      .join(',');

    csvString += '\n';
  }

  return csvString;
}

export function downloadAggCaseData({
  groupKey,
  dnaOrAa,
  coordinateMode,
  dataAggGroup,

  ...rest
}) {
  // console.log(groupKey, dnaOrAa, dataAggGroup);

  let csvString = '';

  // Get list of changing positions.
  // We don't need to pass the store variable in. the positions will be
  // encoded onto each row of dataAggGroup, so just grab one row
  // to see what's there
  const changingPositions = [];
  Object.keys(dataAggGroup[0]).forEach((key) => {
    if (key.slice(0, 4) === 'pos_') {
      // Changing positions are 0-indexed
      changingPositions.push(parseInt(key.slice(4)));
    }
  });

  // If we're in lineage mode, then we need to get SNPs for this lineage
  if (Object.keys(config.group_cols).includes(groupKey)) {
    csvString = downloadAggCaseDataGroup({
      groupKey,
      dataAggGroup,
      changingPositions,
      coordinateMode,

      ...rest,
    });
  } else if (groupKey === GROUP_SNV) {
    csvString = downloadAggCaseDataSnp(
      dnaOrAa,
      dataAggGroup,
      changingPositions
    );
  }

  let blob = new Blob([csvString]);
  let url = URL.createObjectURL(blob);

  return {
    blobURL: url,
  };
}
