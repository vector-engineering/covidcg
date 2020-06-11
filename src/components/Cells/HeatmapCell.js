import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { reds } from '../../utils/colors';

const numColors = reds.length;

const HeatmapCellDiv = styled.div`
  width: 100%;
  height: 100%;
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-start;
  padding: 0px 8px;
  font-weight: normal;
`;

const HeatmapCell = ({ value, min, max, percent }) => {
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
    <HeatmapCellDiv style={{ backgroundColor: color }}>{value}</HeatmapCellDiv>
  );
};

HeatmapCell.propTypes = {
  value: PropTypes.number,
  min: PropTypes.number,
  max: PropTypes.number,
  percent: PropTypes.bool,
};

export default HeatmapCell;
