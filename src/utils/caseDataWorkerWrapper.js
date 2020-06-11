// eslint-disable-next-line import/default
import Worker from './caseData.worker.js';

const caseDataWorker = new Worker();

export const processCaseData = (caseData, callback) => {
  caseDataWorker.onmessage = (e) => {
    callback(JSON.parse(e.data));
  };

  caseData.type = 'processCaseData';

  caseDataWorker.postMessage(JSON.stringify(caseData)); // Send data to our worker.
};

export const aggCaseDataByGroup = (
  { caseData, selectedGene, groupKey, dnaOrAa, dateRange },
  callback
) => {
  caseDataWorker.postMessage(
    JSON.stringify({
      type: 'aggCaseDataByGroup',
      caseData,
      selectedGene,
      groupKey,
      dnaOrAa,
      dateRange,
    })
  );

  caseDataWorker.onmessage = (e) => {
    callback(JSON.parse(e.data));
  };
};
