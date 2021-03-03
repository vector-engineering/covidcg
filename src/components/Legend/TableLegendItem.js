import React, { useState, useEffect } from 'react';
import { observer } from 'mobx-react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import _ from 'underscore';
import { GROUP_SNV } from '../../constants/defs.json';
import { useStores } from '../../stores/connect';
import { formatSnv } from '../../utils/snpUtils';
import { reds } from '../../constants/colors';

const numColors = reds.length;

const Container = styled.div`
  border-bottom: 1px solid #eee;
  display: flex;
  align-items: center;
  ${({ hovered, selected }) => {
    if (selected) return 'background-color: rgba(0,0,0,0.3);';
    if (hovered) return 'background-color: rgba(0,0,0,0.1);';
  }}
  font-size: 12px;
`;

const ColorBar = styled.div`
  height: 100%;
  border-right: 5px solid ${({ color }) => color};
  margin-right: 6px;
`;

const GroupNameContainer = styled.div`
  border-right: 1px #eee solid;
  width: 40%;
  height: 100%;
  display: flex;
  align-items: center;

  font-size: 0.85em;
  line-height: normal;
  word-break: break-all;
`;

const PercentageContainer = styled.div`
  width: 30%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: flex-start;

  padding: 0px 5px;
  font-size: 0.7rem;
`;

const PercentageCell = ({ value, min, max, percent }) => {
  // Find the color for this value
  let color = '#FFFFFF';
  // Add a bit extra since sometimes rounding errors can exclude the max value
  let interval = (max - min) / numColors + 0.00001;

  for (let i = 0; i < numColors; i++) {
    if (value >= min + i * interval && value <= min + (i + 1) * interval) {
      color = reds[i];
      break;
    }
  }

  // Don't show NaNs
  if (Number.isNaN(value) || value === null) {
    value = '';
    color = 'transparent';
  }
  // Format percentages
  else if (percent === true) {
    value = (value * 100).toFixed(2) + '%';
  }

  return (
    <PercentageContainer style={{ backgroundColor: color }}>
      {value}
    </PercentageContainer>
  );
};

PercentageCell.propTypes = {
  value: PropTypes.number,
  min: PropTypes.number,
  max: PropTypes.number,
  percent: PropTypes.bool,
};

const TableLegendItem = observer(
  ({ item, style, maxCounts, updateHoverGroup, updateSelectGroup }) => {
    const { configStore } = useStores();
    const [hovered, setHovered] = useState();
    const [selected, setSelected] = useState();

    const onMouseMove = () => {
      setHovered(true);
      if (item.group !== configStore.hoverGroup) {
        updateHoverGroup(item.group);
      }
    };

    const onMouseOut = () => {
      setHovered(false);
      if (configStore.hoverGroup !== null) {
        updateHoverGroup(null);
      }
    };

    useEffect(() => {
      const _hovered =
        configStore.hoverGroup === null
          ? false
          : configStore.hoverGroup === item.group;

      if (_hovered !== hovered) {
        // console.log(hovered, state.hovered);
        setHovered(hovered);
      }
    }, [configStore.hoverGroup]);

    useEffect(() => {
      let _selected = null;
      if (configStore.selectedGroups.length > 0) {
        if (
          _.findWhere(configStore.selectedGroups, { group: item.group }) !==
          undefined
        ) {
          _selected = true;
        } else {
          _selected = false;
        }
      }

      setSelected(_selected);
    }, [configStore.selectedGroups]);

    return (
      <Container
        style={style}
        onMouseEnter={onMouseMove}
        onMouseOut={onMouseOut}
        onClick={updateSelectGroup}
        hovered={hovered}
        selected={selected}
      >
        <ColorBar color={item.color} />
        <GroupNameContainer data-group={item.group}>
          {configStore.groupKey === GROUP_SNV
            ? formatSnv(item.group, configStore.dnaOrAa)
            : item.group}
        </GroupNameContainer>
        <PercentageCell
          data-group={item.group}
          min={0}
          max={maxCounts}
          value={item.counts}
          percent={false}
        />
        <PercentageCell
          data-group={item.group}
          min={0}
          max={1}
          value={item.percent}
          percent={true}
        />
      </Container>
    );
  }
);
TableLegendItem.propTypes = {
  item: PropTypes.object.isRequired,
  style: PropTypes.object,
  maxCounts: PropTypes.number,
  updateHoverGroup: PropTypes.func,
  updateSelectGroup: PropTypes.func,
};

export default TableLegendItem;
