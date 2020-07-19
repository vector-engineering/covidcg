import React from 'react';
import PropTypes from 'prop-types';
import { toJS } from 'mobx';
import { observer } from 'mobx-react';
import { useStores } from '../stores/connect';
import _ from 'underscore';

import { Vega } from 'react-vega';

// import areaStackSpecInitial from '../vega/area_stack.vl.json';
import barStackSpecInitial from '../vega/bar_stack_v1.vg.json';

const VegaWrapper = observer(({ width }) => {
  // const caseData = toJS(data.case_data);
  const { covidStore } = useStores();

  const handleBrush = (...args) => {
    let dateRange = args[1];
    if (dateRange !== null) {
      covidStore.selectDateRange([
        dateRange[0].getTime(),
        dateRange[1].getTime(),
      ]);
    } else {
      // Reset time range
      covidStore.selectDateRange([-1, -1]);
    }
  };

  // let newCaseData;

  // if (covidStore.groupsToKeep) {
  //   newCaseData = [];
  //   const addedOthersByDate = {};

  //   covidStore.caseData.forEach((row, key) => {
  //     if (!covidStore.groupsToKeep[row.group]) {
  //       row.group = 'other';
  //       row.color = '#cccccc';
  //       if (addedOthersByDate[row.date]) {
  //         covidStore.caseData[addedOthersByDate[row.date]].cases_sum +=
  //           row.cases_sum;
  //       } else {
  //         addedOthersByDate[row.date] = key;
  //         newCaseData.push(row);
  //       }
  //     } else {
  //       newCaseData.push(row);
  //     }
  //   });
  // }

  // Make a deep copy of the Vega spec so we can edit it
  const barStackSpec = JSON.parse(JSON.stringify(barStackSpecInitial));

  // Set the width manually
  barStackSpec['width'] = width;
  barStackSpec['marks'][0]['encode']['enter']['width']['value'] = width;

  // if (areaStackMode === 'percentages') {
  //   areaStackSpec['vconcat'][0]['encoding']['y']['stack'] = 'normalize';
  // }

  // // Adapt labels to groupings
  // if (covidStore.groupKey === 'lineage') {
  //   // y-axis title
  //   areaStackSpec['vconcat'][0]['encoding']['y']['axis']['title'] =
  //     (areaStackMode === 'percentages' ? 'Percent ' : '') +
  //     'Sequences by Lineage';
  //   // Tooltip title
  //   areaStackSpec['vconcat'][0]['encoding']['tooltip'][0]['title'] = 'Lineage';
  // } else if (covidStore.groupKey === 'snp') {
  //   // y-axis title
  //   areaStackSpec['vconcat'][0]['encoding']['y']['axis']['title'] =
  //     (areaStackMode === 'percentages' ? 'Percent ' : '') +
  //     'Sequences by ' +
  //     (covidStore.dnaOrAa === 'dna' ? 'NT' : 'AA') +
  //     ' SNP';
  //   // Tooltip title
  //   areaStackSpec['vconcat'][0]['encoding']['tooltip'][0]['title'] =
  //     (covidStore.dnaOrAa === 'dna' ? 'NT' : 'AA') + ' SNP';
  // }

  let caseData = JSON.parse(JSON.stringify(covidStore.caseData));

  // Compute counts for each group
  let countsPerGroup = {};
  caseData.forEach((row) => {
    if (!Object.prototype.hasOwnProperty.call(countsPerGroup, row.group)) {
      countsPerGroup[row.group] = 0;
    }
    countsPerGroup[row.group] += row.cases_sum;
  });

  // Map counts per group back onto main dataset
  caseData.forEach((row) => {
    row['group_counts'] = countsPerGroup[row.group];
  });

  // For development in Vega Editor
  // console.log(JSON.stringify(caseData));

  return (
    <Vega
      data={{
        cases_by_date_and_group: caseData,
      }}
      spec={barStackSpec}
      signalListeners={{
        detailDomain: _.debounce(handleBrush, 500),
      }}
    />
  );
});

VegaWrapper.propTypes = {
  width: PropTypes.number.isRequired,
};

export default VegaWrapper;
