import React, { useState, useEffect } from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import _ from 'underscore';

import { useStores } from '../../stores/connect';
import SkeletonElement from '../Common/SkeletonElement';
import { asyncStates } from '../../stores/UIStore';
import { mergeLegendItemsIntoOther } from './utils';
import { lighten, transparentize, meetsContrastGuidelines } from 'polished';

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
      return '0px 4px 0px 1px';
    } else {
      return '1px 5px 1px 2px';
    }
  }};
  background-color: ${({ color, hovered, selected }) => {
    if (hovered) {
      return lighten(0.1, color);
    } else if (selected === null) {
      return color;
    } else if (selected) {
      return color;
    } else {
      return transparentize(0.7, color);
    }
  }};

  font-size: 12px;
  font-weight: 500;
  color: ${({ textColor, hovered, selected }) => {
    if (hovered) {
      return transparentize(0.1, textColor);
    } else if (selected === null) {
      return textColor;
    } else if (selected) {
      return textColor;
    } else {
      // Make the blacks a lot lighter
      return textColor === '#fff'
        ? transparentize(0.5, textColor)
        : transparentize(0.8, textColor);
    }
  }};
`;

LegendItem.defaultProps = {
  hovered: false,
  selected: null,
  color: '#ccc',
  textColor: '#fff',
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

const LegendItemWrapper = observer(({ group, color, onItemSelect }) => {
  const { configStore } = useStores();
  const [state, setState] = useState({
    hovered: false,
    selected: false,
  });

  const updateHoverGroup = (hoverGroup) => {
    // Don't fire the action if there's no change
    if (hoverGroup === configStore.hoverGroup) {
      return;
    }
    configStore.updateHoverGroup(hoverGroup);
  };

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

  const handleItemSelect = (e) => {
    onItemSelect(group, e);
  };

  useEffect(() => {
    const hovered =
      configStore.hoverGroup === null
        ? false
        : configStore.hoverGroup === group;
    setState({ ...state, hovered });
  }, [configStore.hoverGroup]);

  useEffect(() => {
    let selected = null;
    if (configStore.selectedGroups.length > 0) {
      if (
        _.findWhere(configStore.selectedGroups, { group: group }) !== undefined
      ) {
        selected = true;
      } else {
        selected = false;
      }
    }
    setState({ ...state, selected });
  }, [configStore.selectedGroups]);

  const scores = meetsContrastGuidelines(color, '#fff');
  const textColor = scores['AALarge'] ? '#fff' : '#000';

  return (
    <LegendItem
      hovered={state.hovered}
      selected={state.selected}
      onMouseEnter={onItemEnter.bind(this, group)}
      onMouseLeave={onItemLeave}
      onMouseDown={handleItemSelect}
      color={color}
      textColor={textColor}
    >
      {group}
    </LegendItem>
  );
});
LegendItemWrapper.propTypes = {
  group: PropTypes.string.isRequired,
  color: PropTypes.string.isRequired,
  onItemSelect: PropTypes.func.isRequired,
};

const VegaLegend = observer(() => {
  const { dataStore, UIStore, configStore } = useStores();

  const [state, setState] = useState({
    legendItems: [],
  });

  const onItemSelect = (selectedGroup, e) => {
    // console.log('select', selectedGroup, e);
    e.preventDefault();

    let newGroups;

    // If the item is already selected, then deselect it
    if (
      _.findWhere(configStore.selectedGroups, { group: selectedGroup }) !==
      undefined
    ) {
      newGroups = _.reject(
        configStore.selectedGroups,
        (group) => group.group == selectedGroup
      );
    } else {
      // Otherwise, add it
      newGroups = [{ group: selectedGroup }];
      // If shift is pressed, then add it to the existing selected groups
      if (UIStore.isKeyPressed(16)) {
        newGroups = newGroups.concat(configStore.selectedGroups);
      }
    }

    configStore.updateSelectedGroups(newGroups);
  };

  const renderLegendKeys = (groupObjs) => {
    return groupObjs.map((obj) => {
      if (!obj.color) {
        return null;
      }

      return (
        <LegendItemWrapper
          key={`legend-item-${obj.group}`}
          color={obj.color}
          group={obj.group}
          onItemSelect={onItemSelect}
        />
      );
    });
  };

  useEffect(() => {
    // Make own copy of the elements, and sort by group
    let legendItems = mergeLegendItemsIntoOther(
      JSON.parse(JSON.stringify(dataStore.caseDataAggGroup)),
      dataStore.groupsToKeep
    );
    // console.log(legendItems, dataStore.caseDataAggGroup);
    legendItems = _.sortBy(legendItems, (row) => row.group);

    setState({ ...state, legendItems: renderLegendKeys(legendItems) });
  }, [dataStore.caseDataAggGroup, dataStore.groupsToKeep]);

  // console.log('RE-RENDERING LEGEND');

  if (UIStore.caseDataState === asyncStates.STARTED) {
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

  return <LegendList>{state.legendItems}</LegendList>;
});

VegaLegend.displayName = 'VegaLegend';

export default VegaLegend;
