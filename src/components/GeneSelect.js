import React, { Component } from 'react';
import PropTypes from 'prop-types'
import _ from 'underscore';

class GeneSelect extends Component {
  constructor(props) {
    super(props);
    this.handleChange = this.handleChange.bind(this);
  }

  handleChange(event) {
    this.props.onChange(event.target.value);
  }

  render() {
    const { genes, value, startPos, endPos } = this.props; 

    // Create option elements
    let option_elements = [];
    genes.forEach(opt => {
      option_elements.push(
        <option key={opt.value} value={opt.value}>{opt.label}</option>
      );
    });

    return (
      <div className='gene-select'>
        <form onSubmit={this.handleSubmit}>
          <label>
            Select a gene to analyze
            <select value={value} onChange={this.handleChange}>
              {option_elements}
            </select>
          </label>
        </form>
        <div className='position-container'>
          <div className='pos-from'>From: {startPos}</div>
          <div className='pos-to'>To: {endPos}</div>
        </div>
      </div>
    );
  }
}

GeneSelect.propTypes = {
  genes: PropTypes.array.isRequired,
  value: PropTypes.string.isRequired,
  startPos: PropTypes.number.isRequired,
  endPos: PropTypes.number.isRequired,
  onChange: PropTypes.func
};

GeneSelect.defaultProps = {
  onChange: () => {}
};

export default GeneSelect;