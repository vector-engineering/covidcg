// Load acknowledgement data

import _ from 'underscore';
import { asyncDataStoreInstance } from '../stores/rootStore';

export function getAckTextsFromAckIds(ackIds) {
  // For each acknowledgement ID, return a object representing the
  // acknowledgement text
  return _.map(ackIds, (ackId) => {
    if (ackId == -1) {
      return {};
    }

    return asyncDataStoreInstance.data.ack_map[ackId];
  });
}
