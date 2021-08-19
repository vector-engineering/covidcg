import React from 'react';
import PropTypes from 'prop-types';

import {
  MarkLegendContainer,
  LegendTitle,
  MarkItemList,
  MarkItemContainer,
  MarkColorBox,
  MarkLabel,
} from './MarkLegend.styles';

const MarkItem = ({ label, color }) => {
  return (
    <MarkItemContainer>
      <MarkColorBox color={color} />
      <MarkLabel>{label}</MarkLabel>
    </MarkItemContainer>
  );
};
MarkItem.propTypes = {
  label: PropTypes.string.isRequired,
  color: PropTypes.string.isRequired,
};

const MarkLegend = ({ title, itemLabels, itemColors }) => {
  const markItems = [];
  itemLabels.forEach((label, i) => {
    markItems.push(
      <MarkItem key={`mark-${i}`} label={label} color={itemColors[i]} />
    );
  });

  return (
    <MarkLegendContainer>
      <LegendTitle>{title}</LegendTitle>
      <MarkItemList>{markItems}</MarkItemList>
    </MarkLegendContainer>
  );
};
MarkLegend.propTypes = {
  title: PropTypes.string.isRequired,
  itemLabels: PropTypes.arrayOf(PropTypes.string).isRequired,
  itemColors: PropTypes.arrayOf(PropTypes.string).isRequired,
};

export default MarkLegend;
