import React, { useState, memo } from 'react';
import styled from 'styled-components';

import { useStores } from '../stores/connect';
import SkeletonElement from './SkeletonElement';
import { asyncStates } from '../stores/uiStore';

const LegendAndButtonWrapper = styled.div`
  display: flex;
  flex-direction: row;
  align-items: flex-start;
  margin-top: 12px;
  margin-bottom: 6px;
  margin-left: 12px;
  margin-right: 24px;
`;

const LegendContainer = styled.div`
  border: 1px #bbb solid;
  padding: 4px 6px;
  border-radius: 2px;
  max-height: 100px;
  overflow-y: scroll;
  transition: max-height 0.3s cubic-bezier(0.22, 1, 0.36, 1);
  ${({ collapsed }) => collapsed && `max-height 5px;`}
`;

const CollapseButton = styled.button`
  border: none;
  cursor: pointer;
  color: #ccc;
  transition: color 0.2s cubic-bezier(0.22, 1, 0.36, 1);
  background-color: rgba(0, 0, 0, 0);
  outline: none;
  width: 32px;
  font-size: 20px;

  &:hover {
    color: black;
  }
`;

const LegendList = styled.div`
  display: flex;
  flex-direction: row;
  flex-wrap: wrap;
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

const VegaLegend = () => {
  const { covidStore, uiStore } = useStores();
  const [collapsed, setCollapsed] = useState(false);

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
        <SkeletonElement delay={1} height={'50px'} />
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

  return (
    <LegendAndButtonWrapper>
      <CollapseButton
        onClick={() => {
          setCollapsed(!collapsed);
        }}
      >
        {collapsed ? '+' : '-'}
      </CollapseButton>
      <LegendContainer collapsed={collapsed}>
        <LegendList>{renderLegendKeys(covidStore.caseDataAggGroup)}</LegendList>
      </LegendContainer>
    </LegendAndButtonWrapper>
  );
};

VegaLegend.displayName = 'VegaLegend';

export default VegaLegend;
