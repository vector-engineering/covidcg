// eslint-disable-next-line import/default
import Worker from './snpData.worker.js';
const snpDataWorker = new Worker();

export const processSelectedSnvs = (data, callback) => {
  snpDataWorker.onmessage = (e) => {
    callback(JSON.parse(e.data));
  };
  data.type = 'processSelectedSnvs';
  snpDataWorker.postMessage(JSON.stringify(data));
};
