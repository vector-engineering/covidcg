// eslint-disable-next-line import/default
import Worker from './mutationData.worker.js';
const snpDataWorker = new Worker();

const callbacks = {};
snpDataWorker.onmessage = (e) => {
  callbacks[e.data.type](e.data);
};

export const processSelectedMutations = (pkg, callback) => {
  const type = 'processSelectedMutations';
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
