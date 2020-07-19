import React from 'react';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import _ from 'underscore';

import { useStores } from '../../stores/connect';
import SkeletonElement from '../SkeletonElement';
import { asyncStates } from '../../stores/uiStore';

const LegendList = styled.div`
  display: flex;
  flex-direction: row;
  flex-wrap: wrap;
  border: 1px #eaeaea solid;
  border-radius: 2px;
  padding: 6px 12px;
`;

const LegendItem = styled.div`
  display: flex;
  align-items: center;
  margin-right: 12px;
`;

const ColorCircle = styled.div`
  ${({ color }) => color && `background-color: ${color};`}
  width: 12px;
  height: 12px;
  border-radius: 50%;
  margin-right: 4px;
`;

const LegendText = styled.span`
  font-size: 10px;
  font-weight: 400;
`;

const VegaLegend = observer(() => {
  const { covidStore, uiStore } = useStores();

  if (uiStore.caseDataState === asyncStates.STARTED) {
    return (
      <div
        style={{
          paddingTop: '12px',
          paddingRight: '24px',
          paddingLeft: '12px',
          paddingBottom: '0px',
        }}
      >
        <SkeletonElement delay={1} height={50} />
      </div>
    );
  }

  const renderLegendKeys = (groupObjs) => {
    return groupObjs.map((obj) => {
      if (!obj.color) {
        return null;
      }
      return (
        <LegendItem key={`${Math.random()}${obj.color}`}>
          <ColorCircle color={obj.color} />
          <LegendText>{obj.group}</LegendText>
        </LegendItem>
      );
    });
  };

  // Make own copy of the elements, and sort by group
  let legendItems = JSON.parse(JSON.stringify(covidStore.caseDataAggGroup));
  legendItems = _.sortBy(legendItems, (row) => row.group);

  return <LegendList>{renderLegendKeys(legendItems)}</LegendList>;
});

VegaLegend.displayName = 'VegaLegend';

export default VegaLegend;
