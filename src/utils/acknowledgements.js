// Load acknowledgement data

import ackIdToAckText from '../../data/acknowledgement_map.json';
import accessionIdToAckId from '../../data/taxon_acknowledgements.json';

import _ from 'underscore';

export function getAckIdsFromAccessionIds(accessionIds) {
  // Do a hash lookup for each accession ID
  return _.map(accessionIds, (row) => {
    return accessionIdToAckId[row['gisaid_id']];
  });
}

export function getAckTextsFromAckIds(ackIds) {
  // For each acknowledgement ID, return a object representing the
  // acknowledgement text
  return _.map(ackIds, (ackId) => {
    return ackIdToAckText[ackId];
  });
}

// accessionIds -> ackTexts
export function getAckTextsFromAccessionIds(accessionIds) {
  return getAckTextsFromAckIds(getAckIdsFromAccessionIds(accessionIds));
}
