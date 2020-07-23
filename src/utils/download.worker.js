import { getAckTextsFromAckIds } from './acknowledgements';
import {
  getDnaSnpsFromGroup,
  getGeneAaSnpsFromGroup,
  getProteinAaSnpsFromGroup,
} from './lineageData';
import { intToISO } from './date';
import _ from 'underscore';

function downloadAcknowledgements({ selectedRows }) {
  let ackIds = _.pluck(selectedRows, 'ack_id');

  // Get the list of selected Accession IDs, and map to
  // acknowledgement texts
  let ackTexts = getAckTextsFromAckIds(ackIds);

  // Write to a CSV string
  // Accession ID and sample date first
  // Then acknowledgement texts
  let csvString =
    'Accession ID,Collection Date,Originating lab,Submitting lab,Authors\n';

  for (let i = 0; i < selectedRows.length; i++) {
    // Write Accession ID
    csvString += selectedRows[i]['Accession ID'] + ',';
    // Write Sample Date
    // Get the date in ISO format, and chop off the time/timezone info at the end
    // So that we get YYYY-MM-DD (the same as the original input format)
    csvString += intToISO(selectedRows[i]['collection_date']) + ',';

    // Write Acknowledgement texts
    // Since these can contain commas, wrap each in double quotes
    csvString += '"' + ackTexts[i]['Originating lab'] + '",';
    csvString += '"' + ackTexts[i]['Submitting lab'] + '",';
    csvString += '"' + ackTexts[i]['Authors'] + '"\n';
  }

  let blob = new Blob([csvString]);
  let url = URL.createObjectURL(blob);

  return {
    blobURL: url,
  };
}

function downloadAggCaseData({
  groupKey,
  dnaOrAa,
  coordinateMode,
  caseDataAggGroup,
}) {
  // console.log(groupKey, dnaOrAa, caseDataAggGroup);

  let csvString = '';

  // Get list of changing positions.
  // We don't need to pass the store variable in. the positions will be
  // encoded onto each row of caseDataAggGroup, so just grab one row
  // to see what's there
  const changingPositions = [];
  Object.keys(caseDataAggGroup[0]).forEach((key) => {
    if (key.slice(0, 4) === 'pos_') {
      // Changing positions are 0-indexed
      changingPositions.push(parseInt(key.slice(4)));
    }
  });

  // If we're in lineage mode, then we need to get SNPs for this lineage
  if (groupKey === 'lineage' || groupKey === 'clade') {
    csvString = downloadAggCaseDataGroup({
      groupKey,
      caseDataAggGroup,
      changingPositions,
      coordinateMode,
    });
  } else if (groupKey === 'snp') {
    csvString = downloadAggCaseDataSnp(
      dnaOrAa,
      caseDataAggGroup,
      changingPositions
    );
  }

  let blob = new Blob([csvString]);
  let url = URL.createObjectURL(blob);

  return {
    blobURL: url,
  };
}

function downloadAggCaseDataGroup({
  groupKey,
  caseDataAggGroup,
  changingPositions,
  coordinateMode,
}) {
  let csvString = '';

  // Write headers
  csvString = 'lineage,seqs,seqs_percent,nt_snps,aa_snps,';
  // Add position column headers
  csvString += _.map(changingPositions, (pos) => (pos + 1).toString()).join(
    ','
  );
  csvString += '\n';

  for (let i = 0; i < caseDataAggGroup.length; i++) {
    let row = caseDataAggGroup[i];
    // Skip if it's the reference row
    if (row['group'] === 'Reference') {
      continue;
    }

    // Write lineage and counts
    csvString +=
      row['group'] + ',' + row['cases_sum'] + ',' + row['cases_percent'] + ',';

    // Get NT SNPs
    let ntSnps = getDnaSnpsFromGroup(groupKey, row['group']);
    // Skip if it's empty
    if (ntSnps.length === 0) {
      csvString += ',';
    } else {
      // Loop thru SNPs and print them as a list
      ntSnps = _.map(ntSnps, (snp) => {
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
    if (coordinateMode === 'gene') {
      aaSnps = getGeneAaSnpsFromGroup(groupKey, row['group']);
    } else if (coordinateMode === 'protein') {
      aaSnps = getProteinAaSnpsFromGroup(groupKey, row['group']);
    }
    // Skip if it's empty
    if (aaSnps.length === 0) {
      csvString += ',';
    } else {
      // Loop thru SNPs and print them as a list
      aaSnps = _.map(aaSnps, (snp) => {
        // Format as gene|pos|ref|alt
        // Position is 0-indexed, so make it 1-indexed
        let label = '';
        if (coordinateMode === 'gene') {
          label = snp['gene'];
        } else if (coordinateMode === 'protein') {
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
    csvString += _.map(
      changingPositions,
      (pos) => row['pos_' + pos.toString()]
    ).join(',');

    csvString += '\n';
  }

  return csvString;
}

function downloadAggCaseDataSnp(dnaOrAa, caseDataAggGroup, changingPositions) {
  let csvString = '';

  // Write headers
  csvString = 'snp_str,';
  // Add DNA headers
  if (dnaOrAa === 'dna') {
    csvString += 'pos,ref,alt,';
  }
  // Add AA headers - just DNA but with gene too
  else {
    csvString += 'gene,pos,ref,alt,';
  }
  // Add count headers
  csvString += 'seqs,seqs_percent,';
  // Add position column headers
  csvString += _.map(changingPositions, (pos) => (pos + 1).toString()).join(
    ','
  );
  csvString += '\n';

  for (let i = 0; i < caseDataAggGroup.length; i++) {
    let row = caseDataAggGroup[i];

    // Write the SNP string
    // For DNA, its pos|ref|alt
    // For AA, its gene|pos|ref|alt
    // And then write the SNP chunks
    if (dnaOrAa === 'dna') {
      // Handle reference row
      if (row['group'] === 'Reference') {
        csvString += 'Reference,,,,';
      } else {
        csvString += [row['pos'], row['ref'], row['alt']].join('|') + ',';
        csvString += [row['pos'], row['ref'], row['alt']].join(',') + ',';
      }
    } else {
      // Handle reference row
      if (row['group'] === 'Reference') {
        csvString += 'Reference,,,,,';
      } else {
        csvString +=
          [row['gene'], row['pos'], row['ref'], row['alt']].join('|') + ',';
        csvString +=
          [row['gene'], row['pos'], row['ref'], row['alt']].join(',') + ',';
      }
    }

    // Write the sequence counts/percents
    csvString += row['cases_sum'] + ',' + row['cases_percent'] + ',';

    // Add letters at positions
    csvString += _.map(
      changingPositions,
      (pos) => row['pos_' + pos.toString()]
    ).join(',');

    csvString += '\n';
  }

  return csvString;
}

self.addEventListener(
  'message',
  function (e) {
    const data = JSON.parse(e.data);
    //console.log('in downloadworker event listener', data);

    let result;
    if (data.type === 'downloadAcknowledgements') {
      // This is a terminal endpoint, we don't need to post a message back
      result = downloadAcknowledgements(data);
    } else if (data.type === 'downloadAggCaseData') {
      result = downloadAggCaseData(data);
    }
    // console.log(result);
    self.postMessage(JSON.stringify(result));
  },
  false
);
