// eslint-disable-next-line import/default
import Worker from './mutationData.worker.js';
const mutationDataWorker = new Worker();

const callbacks = {};
mutationDataWorker.onmessage = (e) => {
  callbacks[e.data.type](e.data);
};

export const processSelectedMutations = (pkg, callback) => {
  const type = 'processSelectedMutations';
  callbacks[type] = callback;
  pkg.type = type;
  mutationDataWorker.postMessage(pkg);
};

export const processCooccurrenceData = (pkg, callback) => {
  const type = 'processCooccurrenceData';
  callbacks[type] = callback;
  pkg.type = type;
  mutationDataWorker.postMessage(pkg);
};
