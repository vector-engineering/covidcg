import React from 'react';
import PropTypes from 'prop-types';

import { StyledLegendList } from './Legend.styles';
import LegendItemWrapper from './LegendItem';

const VegaLegend = ({ legendItems, updateHoverGroup, updateSelectGroup }) => {
  const renderLegendKeys = () => {
    return legendItems.map((obj) => {
      if (!obj.color) {
        return null;
      }

      return (
        <LegendItemWrapper
          key={`legend-item-${obj.group}`}
          color={obj.color}
          group={obj.group}
          updateHoverGroup={updateHoverGroup}
        />
      );
    });
  };

  return (
    <StyledLegendList onMouseDown={updateSelectGroup}>
      {renderLegendKeys()}
    </StyledLegendList>
  );
};

VegaLegend.displayName = 'VegaLegend';
VegaLegend.propTypes = {
  legendItems: PropTypes.array,
  updateHoverGroup: PropTypes.func,
  updateSelectGroup: PropTypes.func,
};

export default VegaLegend;
