import Worker from './caseData.worker.js';

const caseDataWorker = new Worker();

export const processCaseData = (caseData, callback) => {
  caseDataWorker.onmessage = (e) => {
    console.log('out of casedata from process', JSON.parse(e.data));
    callback(JSON.parse(e.data));
  };

  console.log('posting: ', caseData);
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
    console.log('out of casedata from agg', JSON.parse(e.data));
    callback(JSON.parse(e.data));
  };
};
