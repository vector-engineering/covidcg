import React, { useState } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
// import { toJS } from 'mobx';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import _ from 'underscore';

import VegaEmbed from '../../react_vega/VegaEmbed';

// import areaStackSpecInitial from '../vega/area_stack.vl.json';
import barStackSpecInitial from '../../vega/bar_stack_v1.vg.json';

const PlotOptions = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-start;
  margin-bottom: 10px;

  .area-stack-title {
    font-size: 1.25em;
    margin-right: 10px;
    padding-right: 10px;
    padding-left: 18px;

    border-right: 1px solid #ccc;
  }
`;

const AreaStackSelectContainer = styled.div`
  font-weight: normal;
  select {
    margin-left: 0.65em;
    padding: 1px 4px;
    border-radius: 3px;
  }
`;

const AreaStackModeSelect = ({ mode, onChange }) => {
  return (
    <AreaStackSelectContainer>
      <label>
        Display mode:
        <select value={mode} onChange={onChange}>
          <option value="counts">Counts</option>
          <option value="percentages">Percentages</option>
        </select>
      </label>
    </AreaStackSelectContainer>
  );
};
AreaStackModeSelect.propTypes = {
  mode: PropTypes.string.isRequired,
  onChange: PropTypes.func.isRequired,
};

const VegaStackedBars = observer(({ width }) => {
  const { covidStore } = useStores();

  // 'percentages' or 'counts'
  const [areaStackMode, setAreaStackMode] = useState('counts');
  const onChangeAreaStackMode = (event) => setAreaStackMode(event.target.value);

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

  const handleHoverGroup = (...args) => {
    // Don't fire the action if there's no change
    let hoverGroup = args[1] === null ? null : args[1]['group'];
    if (hoverGroup === covidStore.hoverGroup) {
      return;
    }
    covidStore.updateHoverGroup(hoverGroup);
  };

  const handleSelected = (...args) => {
    // console.log(args);
    covidStore.updateSelectedGroups(args[1]);
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

  if (areaStackMode === 'percentages') {
    let caseDataSpec = _.findWhere(barStackSpec['data'], {
      name: 'cases_by_date_and_group',
    });
    let stackTransform = _.findWhere(caseDataSpec['transform'], {
      type: 'stack',
    });
    stackTransform['offset'] = 'normalize';
  }

  // Adapt labels to groupings
  let detailGroup = _.findWhere(barStackSpec['marks'], { name: 'detail' });
  // let detailBars = detailGroup['marks'][0]['marks'][0];
  if (covidStore.groupKey === 'lineage') {
    // y-axis title
    detailGroup['axes'][1]['title'] =
      (areaStackMode === 'percentages' ? 'Percent ' : '') +
      'Sequences by Lineage';
  } else if (covidStore.groupKey === 'snp') {
    // y-axis title
    detailGroup['axes'][1]['title'] =
      (areaStackMode === 'percentages' ? 'Percent ' : '') +
      'Sequences by ' +
      (covidStore.dnaOrAa === 'dna' ? 'NT' : 'AA') +
      ' SNP';
  }

  let caseData = JSON.parse(JSON.stringify(covidStore.caseData));
  let selectedGroups = JSON.parse(JSON.stringify(covidStore.selectedGroups));

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

  let areaStackTitle = 'Lineage ';
  if (covidStore.groupKey === 'lineage') {
    areaStackTitle = 'Lineage ';
  } else if (covidStore.groupKey === 'snp') {
    if (covidStore.dnaOrAa === 'dna') {
      areaStackTitle = 'NT';
    } else {
      areaStackTitle = 'AA';
    }
    areaStackTitle += ' SNP ';
  }
  areaStackTitle += areaStackMode === 'percentages' ? 'Percentages' : 'Counts';
  areaStackTitle += ' Over Time';

  return (
    <div>
      <PlotOptions>
        <span className="area-stack-title">{areaStackTitle}</span>
        <AreaStackModeSelect
          mode={areaStackMode}
          onChange={onChangeAreaStackMode}
        />
      </PlotOptions>

      <div style={{ width: `${width - 150}px` }}>
        <VegaEmbed
          data={{
            cases_by_date_and_group: caseData,
            selected: selectedGroups,
          }}
          spec={barStackSpec}
          signalListeners={{
            detailDomain: _.debounce(handleBrush, 500),
            hoverBar: _.throttle(handleHoverGroup, 100),
          }}
          dataListeners={{
            selected: handleSelected,
          }}
          signals={{
            hoverBar: covidStore.hoverGroup,
          }}
          width={width}
        />
      </div>
    </div>
  );
});

VegaStackedBars.propTypes = {
  width: PropTypes.number.isRequired,
};

export default VegaStackedBars;
