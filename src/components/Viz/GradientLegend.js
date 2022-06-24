import React from 'react';
import PropTypes from 'prop-types';

import {
  LegendContainer,
  LegendTitle,
  LegendGradientBar,
  LegendAxis,
  TickContainer,
  TickLine,
  TickLabel
} from './GradientLegend.styles';

const TickItem = ({ loc, label }) => {
  return (
    <TickContainer loc={loc}>
      <TickLine />
      <TickLabel>{label}</TickLabel>
    </TickContainer>
  );
};
TickItem.propTypes = {
  loc: PropTypes.number.isRequired,
  label: PropTypes.string.isRequired
};

const GradientLegend = ({ title, gradient, ticks, tickLabels }) => {

  const tickItems = [];
  ticks.forEach((tick, i) => {
    tickItems.push(
      <TickItem key={`tick-${i}`} loc={tick} label={tickLabels[i]} />
    );
  });

  return (
    <LegendContainer>
      <LegendTitle>{title}</LegendTitle>
      <LegendGradientBar gradient={gradient} />
      <LegendAxis>{tickItems}</LegendAxis>
    </LegendContainer>
  );
};
GradientLegend.propTypes = {
  title: PropTypes.string.isRequired,
  gradient: PropTypes.string.isRequired,
  ticks: PropTypes.arrayOf(PropTypes.number),
  tickLabels: PropTypes.arrayOf(PropTypes.string),
};

export default GradientLegend;
