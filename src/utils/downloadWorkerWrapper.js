// eslint-disable-next-line import/default
import Worker from './download.worker.js';
const downloadWorker = new Worker();

export const downloadAcknowledgementsData = (data, callback) => {
  downloadWorker.onmessage = (e) => {
    callback(JSON.parse(e.data));
  };
  data.type = 'downloadAcknowledgementsData';
  downloadWorker.postMessage(JSON.stringify(data));
};

export const downloadAggCaseData = (data, callback) => {
  downloadWorker.onmessage = (e) => {
    callback(JSON.parse(e.data));
  };
  data.type = 'downloadAggCaseData';
  downloadWorker.postMessage(JSON.stringify(data));
};

export const downloadAccessionIdsData = (data, callback) => {
  downloadWorker.onmessage = (e) => {
    callback(JSON.parse(e.data));
  };
  data.type = 'downloadAccessionIdsData';
  downloadWorker.postMessage(JSON.stringify(data));
};
