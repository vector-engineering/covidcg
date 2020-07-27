import React, { useState, useEffect } from 'react';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import _ from 'underscore';

import { useStores } from '../../stores/connect';
import SkeletonElement from '../Common/SkeletonElement';
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

  border: ${({ hovered, selected }) => {
    if (hovered) {
      return '1px solid #666';
    } else if (selected) {
      return '1px solid #000';
    } else {
      return 'none';
    }
  }};
  border-radius: 3px;
  margin: ${({ hovered, selected }) => {
    if (hovered || selected) {
      return '0px 0px 0px 0px';
    } else {
      return '1px 1px 1px 1px';
    }
  }};
  background-color: ${({ hovered, selected }) => {
    if (hovered) {
      return 'rgba(0,0,0,0.05)';
    } else if (selected) {
      return 'rgba(0,0,0,0.1)';
    } else {
      return 'transparent';
    }
  }};
`;

LegendItem.defaultProps = {
  hovered: false,
  selected: null,
};

const ColorCircle = styled.div`
  background-color: ${({ color, selected }) => {
    if (selected === null || selected) {
      return color;
    } else {
      return '#CCC';
    }
  }};
  width: 12px;
  height: 12px;
  border-radius: 50%;
  margin-right: 4px;
`;

ColorCircle.defaultProps = {
  color: '#000',
  selected: null,
};

const LegendText = styled.span`
  font-size: 10px;
  font-weight: 400;
`;

const VegaLegend = observer(() => {
  const { covidStore, uiStore } = useStores();
  const [shiftKeyPressed, setShiftKeyPressed] = useState(false);

  const onKeyDown = (e) => {
    // Shift key = 16
    if (e.keyCode === 16) {
      setShiftKeyPressed(true);
    }
  };

  const onKeyUp = (e) => {
    if (e.keyCode === 16) {
      setShiftKeyPressed(false);
    }
  };

  useEffect(() => {
    document.addEventListener('keydown', onKeyDown, false);
    document.addEventListener('keyup', onKeyUp, false);

    return () => {
      document.removeEventListener('keydown', onKeyDown, false);
      document.removeEventListener('keyup', onKeyUp, false);
    };
  });

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

  const onItemSelect = (selectedGroup, e) => {
    // console.log('select', selectedGroup, e);
    e.preventDefault();

    let newGroups;

    // If the item is already selected, then deselect it
    if (
      _.findWhere(covidStore.selectedGroups, { group: selectedGroup }) !==
      undefined
    ) {
      newGroups = _.reject(
        covidStore.selectedGroups,
        (group) => group.group == selectedGroup
      );
    } else {
      // Otherwise, add it
      newGroups = [{ group: selectedGroup }];
      // If shift is pressed, then add it to the existing selected groups
      if (shiftKeyPressed) {
        newGroups = newGroups.concat(covidStore.selectedGroups);
      }
    }

    covidStore.updateSelectedGroups(newGroups);
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

      let itemSelected = null;
      if (covidStore.selectedGroups.length > 0) {
        if (
          _.findWhere(covidStore.selectedGroups, { group: obj.group }) !==
          undefined
        ) {
          itemSelected = true;
        } else {
          itemSelected = false;
        }
      }

      return (
        <LegendItem
          key={`${Math.random()}${obj.color}`}
          hovered={covidStore.hoverGroup === obj.group}
          selected={itemSelected}
          onMouseEnter={onItemEnter.bind(this, obj.group)}
          onMouseLeave={onItemLeave}
          onMouseDown={onItemSelect.bind(this, obj.group)}
        >
          <ColorCircle color={obj.color} selected={itemSelected} />
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
