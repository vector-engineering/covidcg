// Download utilities

import { getLocationNameByIds } from './location';
import { intToISO } from './date';

export function downloadBlobURL(blob_url, filename) {
  const link = window.document.getElementById('hidden-download-link');
  link.setAttribute('href', blob_url);
  link.setAttribute('download', filename);
  link.click();
}

// Generate a string, for downloaded file names, that represents
// the current selection of sequences
export function generateSelectionString(
  prefix,
  suffix,
  groupKey,
  dnaOrAa,
  selectedGene,
  selectedLocationIds,
  dateRange
) {
  // Prepare some regex stuff
  const whitespace_pattern = /\s+/g;

  let out = prefix + '_';

  // Add group key
  out += groupKey + '_';
  // Add DNA/AA, only if in SNP mode
  if (groupKey === 'snp') {
    out += dnaOrAa + '_';
  }
  // Add selected gene
  out += selectedGene.gene.replace(whitespace_pattern, '') + '_';

  // Add locations
  let locNames = getLocationNameByIds(selectedLocationIds);
  console.log(locNames);
  // Don't add them as text, that will take up way too much space
  // Need to find a better way to embed location here

  // Add date range
  let startDate =
    dateRange[0] === -1 ? new Date('2020-01-01').getTime() : dateRange[0];
  let endDate =
    dateRange[1] === -1 ? new Date('2020-12-31').getTime() : dateRange[1];

  out += intToISO(startDate) + '-';
  out += intToISO(endDate) + '_';

  out += '.' + suffix;

  return out;
}
