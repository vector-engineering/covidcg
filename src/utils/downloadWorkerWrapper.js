// eslint-disable-next-line import/default
import Worker from './download.worker.js';
const downloadWorker = new Worker();

export const downloadAcknowledgements = (data, callback) => {
  downloadWorker.onmessage = (e) => {
    callback(JSON.parse(e.data));
  };

  //console.log('Download worker posting:', data);
  data.type = 'downloadAcknowledgements';

  downloadWorker.postMessage(JSON.stringify(data));
};
