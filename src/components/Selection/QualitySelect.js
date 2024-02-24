import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { config } from '../../config';

import QuestionButton from '../Buttons/QuestionButton';

import {
  QualitySelectContainer,
  FormRow,
  TitleColumn,
  FormColumn,
} from './QualitySelect.styles';

const QualitySelect = observer(
  ({ sequenceLengthRange, percentAmbiguousRange, updateQualityFilters }) => {
    // Only render this if we have the quality filters available
    if (config['virus'] !== 'sars2') {
      return '';
    }

    const handleChange = (field, position, event) => {
      //console.log(field, position, event.target.value);
      let rng;
      if (field === 'sequenceLengthRange') {
        rng = sequenceLengthRange;
      } else if (field === 'percentAmbiguousRange') {
        rng = percentAmbiguousRange;
      }
      rng[position] =
        event.target.value === '' ? null : parseFloat(event.target.value);

      updateQualityFilters(field, rng);
    };

    const qualityFilterItems = [];

    qualityFilterItems.push(
      <FormRow key={`quality-filter-length`}>
        <TitleColumn>Sequence Length (bases)</TitleColumn>
        <FormColumn>
          <label htmlFor={`quality-filter-length-start`}>Minimum</label>
          <input
            type="number"
            id={`quality-filter-length-start`}
            name="quality-filter-length-start"
            value={
              sequenceLengthRange[0] === null ? '' : sequenceLengthRange[0]
            }
            min={0}
            step={1}
            onChange={handleChange.bind(this, 'sequenceLengthRange', 0)}
          />
        </FormColumn>
        <FormColumn>
          <label htmlFor={`quality-filter-length-end`}>Maximum</label>
          <input
            type="number"
            id={`quality-filter-length-end`}
            name="quality-filter-length-end"
            value={
              sequenceLengthRange[1] === null ? '' : sequenceLengthRange[1]
            }
            min={0}
            step={1}
            onChange={handleChange.bind(this, 'sequenceLengthRange', 1)}
          />
        </FormColumn>
      </FormRow>
    );

    qualityFilterItems.push(
      <FormRow key={`quality-filter-percent-ambiguous`}>
        <TitleColumn>% Ambiguous (% N)</TitleColumn>
        <FormColumn>
          <label htmlFor={`quality-filter-percent-ambiguous-start`}>
            Minimum
          </label>
          <input
            type="number"
            id={`quality-filter-percent-ambiguous-start`}
            name="quality-filter-percent-ambiguous-start"
            value={
              percentAmbiguousRange[0] === null ? '' : percentAmbiguousRange[0]
            }
            min={0}
            onChange={handleChange.bind(this, 'percentAmbiguousRange', 0)}
          />
        </FormColumn>
        <FormColumn>
          <label htmlFor={`quality-filter-percent-ambiguous-end`}>
            Maximum
          </label>
          <input
            type="number"
            id={`quality-filter-percent-ambiguous-end`}
            name="quality-filter-percent-ambiguous-end"
            value={
              percentAmbiguousRange[1] === null ? '' : percentAmbiguousRange[1]
            }
            min={0}
            onChange={handleChange.bind(this, 'percentAmbiguousRange', 1)}
          />
        </FormColumn>
      </FormRow>
    );

    return (
      <QualitySelectContainer>
        {' '}
        <span className="title">
          Sequence Quality
          <QuestionButton
            data-tip="<p>Filter sequences based on their quality</p>"
            data-html={true}
            data-place="left"
            data-for="main-tooltip"
          />
        </span>
        {qualityFilterItems}
      </QualitySelectContainer>
    );
  }
);

QualitySelect.propTypes = {
  sequenceLengthRange: PropTypes.arrayOf(PropTypes.number),
  percentAmbiguousRange: PropTypes.arrayOf(PropTypes.number),
  updateQualityFilters: PropTypes.func,
};

export default QualitySelect;
