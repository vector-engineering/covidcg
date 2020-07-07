import React from 'react';
import { VegaLite } from 'react-vega';
import { toJS } from 'mobx';
import { useStores } from '../stores/connect';

const VegaWrapper = ({ data, spec, signalListeners }) => {
  const caseData = toJS(data.case_data);
  const { covidStore } = useStores();

  let newCaseData;

  if (covidStore.groupsToKeep) {
    newCaseData = [];
    const addedOthersByDate = {};

    caseData.forEach((row, key) => {
      if (!covidStore.groupsToKeep[row.group]) {
        row.group = 'other';
        row.color = '#cccccc';
        if (addedOthersByDate[row.date]) {
          caseData[addedOthersByDate[row.date]].cases_sum += row.cases_sum;
        } else {
          addedOthersByDate[row.date] = key;
          newCaseData.push(row);
        }
      } else {
        newCaseData.push(row);
      }
    });
  }

  return (
    <VegaLite
      data={{ case_data: newCaseData || data.case_data }}
      spec={spec}
      signalListeners={signalListeners}
    />
  );
};

export default VegaWrapper;
