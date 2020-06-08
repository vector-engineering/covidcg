import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';

const reds = [
  '#FFF5F0',
  '#FEF1EB',
  '#FEEEE6',
  '#FEEAE1',
  '#FEE7DC',
  '#FEE3D7',
  '#FEE0D2',
  '#FDDACB',
  '#FDD4C3',
  '#FDCEBB',
  '#FCC8B3',
  '#FCC2AB',
  '#FCBCA3',
  '#FCB59B',
  '#FCAF93',
  '#FCA88B',
  '#FCA184',
  '#FC9B7C',
  '#FC9474',
  '#FB8D6D',
  '#FB8767',
];

const numColors = reds.length;

const HeatmapCellDiv = styled.div`
  width: 100%;
  height: 100%;
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: center;
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
  if (Number.isNaN(value)) {
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
