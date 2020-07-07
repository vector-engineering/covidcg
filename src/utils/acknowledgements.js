// Load acknowledgement data

import ackIdToAckText from '../../data/acknowledgement_map.json';

import _ from 'underscore';

export function getAckTextsFromAckIds(ackIds) {
  // For each acknowledgement ID, return a object representing the
  // acknowledgement text
  return _.map(ackIds, (ackId) => {
    return ackIdToAckText[ackId];
  });
}
