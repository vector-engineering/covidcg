import React, { useState, useEffect } from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { meetsContrastGuidelines } from 'polished';
import _ from 'underscore';

import { GROUP_KEYS } from '../../constants/config';
import { formatSnv } from '../../utils/snpUtils';
import { LegendItem } from './Legend.styles';

const LegendItemWrapper = observer(({ group, color, updateHoverGroup }) => {
  const { configStore } = useStores();
  const [state, setState] = useState({
    text:
      configStore.groupKey === GROUP_KEYS.GROUP_SNV
        ? formatSnv(group, configStore.dnaOrAa)
        : group,
    hovered: false,
    selected: false,
    textColor: meetsContrastGuidelines(color, '#fff')['AALarge']
      ? '#fff'
      : '#000',
  });

  const onMouseMove = () => {
    setState({
      ...state,
      hovered: true,
    });
    if (group !== configStore.hoverGroup) {
      updateHoverGroup(group);
    }
  };

  const onMouseOut = () => {
    setState({
      ...state,
      hovered: false,
    });
    if (configStore.hoverGroup !== null) {
      updateHoverGroup(null);
    }
  };

  useEffect(() => {
    const hovered =
      configStore.hoverGroup === null
        ? false
        : configStore.hoverGroup === group;

    if (hovered !== state.hovered) {
      // console.log(hovered, state.hovered);
      setState({ ...state, hovered });
    }
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

  // console.log('re-rendering legend item');

  return (
    <LegendItem
      hovered={state.hovered}
      selected={state.selected}
      color={color}
      textColor={state.textColor}
      data-group={group}
      onMouseEnter={onMouseMove}
      onMouseOut={onMouseOut}
    >
      {state.text}
    </LegendItem>
  );
});
LegendItemWrapper.propTypes = {
  group: PropTypes.string.isRequired,
  color: PropTypes.string.isRequired,
};

export default LegendItemWrapper;
