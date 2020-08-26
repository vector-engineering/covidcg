// eslint-disable-next-line import/default
import Worker from './snpData.worker.js';
const snpDataWorker = new Worker();

const callbacks = {};
snpDataWorker.onmessage = (e) => {
  callbacks[e.data.type](e.data);
};

export const processSelectedSnvs = (pkg, callback) => {
  const type = 'processSelectedSnvs';
  callbacks[type] = callback;
  pkg.type = type;
  snpDataWorker.postMessage(pkg);
};

export const processCooccurrenceData = (pkg, callback) => {
  const type = 'processCooccurrenceData';
  callbacks[type] = callback;
  pkg.type = type;
  snpDataWorker.postMessage(pkg);
};
