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
  padding: 0px 6px;

  border: ${(props) => (props.hovered ? '1px solid #666' : 'none')};
  margin: ${(props) => (props.hovered ? '0px 0px 0px 0px' : '1px 1px 1px 1px')};
`;

LegendItem.defaultProps = {
  hovered: false,
};

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

  const onItemEnter = (hoverGroup, e) => {
    // console.log('enter', hoverGroup, e);
    e.preventDefault();
    updateHoverGroup(hoverGroup);
  };

  const onItemLeave = (e) => {
    // console.log('leave', e);
    e.preventDefault();
    updateHoverGroup(null);
  };

  const updateHoverGroup = (hoverGroup) => {
    // Don't fire the action if there's no change
    if (hoverGroup === covidStore.hoverGroup) {
      return;
    }
    covidStore.updateHoverGroup(hoverGroup);
  };

  const renderLegendKeys = (groupObjs) => {
    return groupObjs.map((obj) => {
      if (!obj.color) {
        return null;
      }
      return (
        <LegendItem
          key={`${Math.random()}${obj.color}`}
          hovered={covidStore.hoverGroup === obj.group}
          onMouseEnter={onItemEnter.bind(this, obj.group)}
          onMouseLeave={onItemLeave}
        >
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
