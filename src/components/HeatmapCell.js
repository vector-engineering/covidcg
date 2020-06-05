import React from 'react';
import PropTypes from 'prop-types';

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

class HeatmapCell extends React.Component {
  constructor(props) {
    super(props);
  }

  render() {
    const { min, max, percent } = this.props;
    let { value } = this.props;

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

    // Format percentages
    if (percent === true) {
      value = (value * 100).toFixed(2) + '%';
    }

    return (
      <div className="heatmap-cell" style={{ backgroundColor: color }}>
        {value}
      </div>
    );
  }
}

HeatmapCell.propTypes = {
  value: PropTypes.number,
  min: PropTypes.number,
  max: PropTypes.number,
  numColors: PropTypes.number,
  percent: PropTypes.bool,
};

export default HeatmapCell;
