import { toJS } from 'mobx';

export const mergeGroupsIntoOther = (data, _groupsToKeep) => {
  let newCaseData;
  const groupsToKeep = toJS(_groupsToKeep);

  if (groupsToKeep) {
    newCaseData = [];
    const addedOthersByDate = {};

    data.forEach((row, key) => {
      if (!groupsToKeep.includes(row.group)) {
        row.group = 'other';
        row.color = '#aaa';
        if (addedOthersByDate[row.date]) {
          data[addedOthersByDate[row.date]].cases_sum += row.cases_sum;
        } else {
          addedOthersByDate[row.date] = key;
          newCaseData.push(row);
        }
      } else {
        newCaseData.push(row);
      }
    });

    return newCaseData;
  }

  throw new Error('mergeGroupsIntoOther: no groups to keep');
};

export const mergeLegendItemsIntoOther = (data, _groupsToKeep) => {
  let newCaseData;
  const groupsToKeep = toJS(_groupsToKeep);

  if (groupsToKeep) {
    newCaseData = [];
    let alreadyPushedOther = false;

    data.forEach((row, key) => {
      if (!groupsToKeep.includes(row.group)) {
        row.group = 'other';
        row.color = '#aaa';

        if (!alreadyPushedOther) {
          newCaseData.push(row);
          alreadyPushedOther = true;
        }
      } else {
        newCaseData.push(row);
      }
    });

    return newCaseData;
  }

  throw new Error('mergeGroupsIntoOther: no groups to keep');
};
