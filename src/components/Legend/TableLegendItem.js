import { observer } from 'mobx-react';
import React, { useState, useEffect } from 'react';
import styled from 'styled-components';
import _ from 'underscore';
import { GROUP_KEYS } from '../../constants/config';
import { useStores } from '../../stores/connect';
import { formatSnv } from '../../utils/snpUtils';

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
  margin-right: 10px;
`;

const TableLegendItem = observer(
  ({ style, group, color, updateHoverGroup, updateSelectGroup }) => {
    const { configStore } = useStores();
    const [hovered, setHovered] = useState();
    const [selected, setSelected] = useState();

    const onMouseMove = () => {
      setHovered(true);
      if (group !== configStore.hoverGroup) {
        updateHoverGroup(group);
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
          : configStore.hoverGroup === group;

      if (_hovered !== hovered) {
        // console.log(hovered, state.hovered);
        setHovered(hovered);
      }
    }, [configStore.hoverGroup]);

    useEffect(() => {
      let _selected = null;
      if (configStore.selectedGroups.length > 0) {
        if (
          _.findWhere(configStore.selectedGroups, { group: group }) !==
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
        onMouseDown={updateSelectGroup}
        data-group={group}
        hovered={hovered}
        selected={selected}
      >
        <ColorBar color={color} />
        {configStore.groupKey === GROUP_KEYS.GROUP_SNV
          ? formatSnv(group, configStore.dnaOrAa)
          : group}
      </Container>
    );
  }
);

export default TableLegendItem;
