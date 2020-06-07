import React, { Component } from 'react';
import PropTypes from 'prop-types';

class GeneSelect extends Component {
  constructor(props) {
    super(props);
    this.handleChange = this.handleChange.bind(this);
  }

  handleChange(event) {
    this.props.onChange(event.target.value);
  }

  render() {
    const { genes, selectedGene } = this.props;

    // Create option elements
    let optionElements = [];
    genes.forEach((opt) => {
      optionElements.push(
        <option key={opt.value} value={opt.value}>
          {opt.label}
        </option>
      );
    });

    return (
      <div className="gene-select">
        <form onSubmit={this.handleSubmit}>
          <label>
            Select a gene to analyze
            <select value={selectedGene.gene} onChange={this.handleChange}>
              {optionElements}
            </select>
          </label>
        </form>
        <div className="position-container">
          <div className="pos-from">From: {selectedGene.start}</div>
          <div className="pos-to">To: {selectedGene.end}</div>
        </div>
      </div>
    );
  }
}

GeneSelect.propTypes = {
  genes: PropTypes.array.isRequired,
  selectedGene: PropTypes.object.isRequired,
  onChange: PropTypes.func,
};

GeneSelect.defaultProps = {
  onChange: () => {},
};

export default GeneSelect;
