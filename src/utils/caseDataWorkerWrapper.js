// eslint-disable-next-line import/default
import Worker from './caseData.worker.js';

const caseDataWorker = new Worker();

export const processCaseData = (pkg, callback) => {
  caseDataWorker.onmessage = (e) => {
    callback(e.data);
  };
  pkg.type = 'processCaseData';
  caseDataWorker.postMessage(pkg); // Send data to our worker.
};

export const aggCaseDataByGroup = (pkg, callback) => {
  caseDataWorker.onmessage = (e) => {
    callback(e.data);
  };
  pkg.type = 'aggCaseDataByGroup';
  caseDataWorker.postMessage(pkg);
};
