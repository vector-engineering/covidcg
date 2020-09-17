// Load acknowledgement data

//import ackIdToAckText from '../../data/ack_map.json';

import _ from 'underscore';

export function getAckTextsFromAckIds(ackIds) {
  // For each acknowledgement ID, return a object representing the
  // acknowledgement text
  return _.map(ackIds, (ackId) => {
    if (ackId == -1) {
      return {};
    }
    return {};
    //return ackIdToAckText[ackId];
  });
}
