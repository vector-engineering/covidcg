import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import _ from 'underscore';

import { Row } from 'react-data-grid';

import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

const RowWrapper = styled.div`
  position: relative;

  .row-cover {
    display: ${({ hovered, selected }) => {
      if (hovered || selected) {
        return 'block';
      } else {
        return 'none';
      }
    }};
    position: absolute;
    pointer-events: none; // allow clicking elements underneath
    top: 0;
    left: 1px;
    width: 100%;
    height: 100%;
    z-index: 3;
    outline: ${({ hovered, selected }) => {
      if (hovered) {
        return '1px solid #666';
      } else if (selected) {
        return '1px solid #000';
      } else {
        return 'none;';
      }
    }};
    background-color: ${({ hovered, selected }) => {
      if (hovered) {
        return 'rgba(0,0,0,0.1)';
      } else if (selected) {
        return 'rgba(0,0,0,0.1)';
      } else {
        return 'none';
      }
    }};
  }

  .rdg-row {
    opacity: ${({ selected }) =>
      selected !== null && !selected ? '0.6' : 'initial'};
  }
`;

RowWrapper.defaultProps = {
  hovered: false,
  selected: null,
};

const RowRenderer = observer(({ row, ...rest }) => {
  const { dataStore } = useStores();
  // console.log(row.group);

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
    if (hoverGroup === dataStore.hoverGroup) {
      return;
    }
    dataStore.updateHoverGroup(hoverGroup);
  };

  let rowSelected = null;
  if (dataStore.selectedGroups.length > 0) {
    if (
      _.findWhere(dataStore.selectedGroups, { group: row.group }) !== undefined
    ) {
      rowSelected = true;
    } else {
      rowSelected = false;
    }
  }

  return (
    <RowWrapper
      hovered={dataStore.hoverGroup === row.group}
      selected={rowSelected}
      onMouseEnter={onItemEnter.bind(this, row.group)}
      onMouseLeave={onItemLeave}
    >
      <div className="row-cover"></div>
      <Row row={row} {...rest} />
    </RowWrapper>
  );
});

RowRenderer.propTypes = {
  row: PropTypes.any,
};

export default RowRenderer;
