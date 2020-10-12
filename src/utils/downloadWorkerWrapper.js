// eslint-disable-next-line import/default
import Worker from './download.worker.js';
const downloadWorker = new Worker();

const callbacks = {};
downloadWorker.onmessage = (e) => {
  callbacks[e.data.type](e.data);
};

// export const downloadAcknowledgementsData = (pkg, callback) => {
//   const type = 'downloadAcknowledgementsData';
//   callbacks[type] = callback;
//   pkg.type = type;
//   downloadWorker.postMessage(pkg);
// };

export const downloadAggCaseData = (pkg, callback) => {
  const type = 'downloadAggCaseData';
  callbacks[type] = callback;
  pkg.type = type;
  downloadWorker.postMessage(pkg);
};

export const downloadAccessionIdsData = (pkg, callback) => {
  const type = 'downloadAccessionIdsData';
  callbacks[type] = callback;
  pkg.type = type;
  downloadWorker.postMessage(pkg);
};
