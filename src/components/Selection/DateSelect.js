import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { intToISO } from '../../utils/date';
import { config } from '../../config';

import {
  Container,
  DateForm,
  FormColumn,
  PresetForm,
} from './DateSelect.styles';

import QuestionButton from '../Buttons/QuestionButton';

const DEFAULT_PRESET_OPTION = 'date-preset-default';
const PRESET_START = 'preset-start';
const presetDateOptions = [
  <option
    key={DEFAULT_PRESET_OPTION}
    value={DEFAULT_PRESET_OPTION}
    disabled={true}
  >
    Select a preset date range
  </option>,
  <option key={PRESET_START} value={PRESET_START}>
    {config.virus === 'sars2'
      ? 'Since pandemic start'
      : 'Date of first sequence'}
  </option>,
  <option key={'separator-1'} disabled={true}>
    ──────────
  </option>,
];

const curPresets = ['This week', 'This month', 'This year'];
curPresets.forEach((presetName) => {
  presetDateOptions.push(
    <option key={presetName} value={presetName}>
      {presetName}
    </option>
  );
});

presetDateOptions.push(
  <option key={'separator-2'} disabled={true}>
    ──────────
  </option>
);

// name: days from today
const lastPresetMap = {
  'Last week': 7,
  'Last 2 weeks': 14,
  'Last month': 30,
  'Last 3 months': 90,
  'Last 6 months': 180,
  'Last year': 365,
};
Object.keys(lastPresetMap).forEach((presetName) => {
  presetDateOptions.push(
    <option key={presetName} value={presetName}>
      {presetName}
    </option>
  );
});

const DateSelect = observer(
  ({ startDate, endDate, updateDateRange, title, hintText }) => {
    const maxDate = new Date().getTime();
    const maxDateISO = intToISO(maxDate);

    const handleStartChange = (event) => {
      updateDateRange(event.target.value, endDate);
    };

    const handleEndChange = (event) => {
      updateDateRange(startDate, event.target.value);
    };

    const handlePresetChange = (event) => {
      const option = event.target.value;

      let startDate = new Date().getTime();
      const endDate = new Date().getTime();

      // If the option chosen was 'Last ...'
      if (Object.keys(lastPresetMap).includes(option)) {
        const daysAgo = lastPresetMap[option];
        // 24 hr/day * 60 min/hr * 60 s/min * 1000 ms/s
        startDate = intToISO(endDate - daysAgo * 24 * 60 * 60 * 1000);
      }
      // If the option chosen was 'This ...'
      else if (curPresets.includes(option)) {
        if (option === 'This week') {
          // Day of the week, starting at Sunday (0)
          const daysAgo = new Date().getDay();
          startDate = intToISO(endDate - daysAgo * 24 * 60 * 60 * 1000);
        } else if (option === 'This month') {
          // Day of the month, starting at 1
          const daysAgo = new Date().getDate() - 1;
          startDate = intToISO(endDate - daysAgo * 24 * 60 * 60 * 1000);
        } else if (option === 'This year') {
          startDate = intToISO(
            new Date(new Date().getFullYear().toString()).getTime()
          );
        }
      } else if (option === PRESET_START) {
        startDate = config.min_date;
      }

      updateDateRange(startDate, intToISO(endDate));
    };

    return (
      <Container>
        <span className="title">
          {title}
          <QuestionButton
            data-tip={hintText}
            data-html={true}
            data-place="bottom"
            data-for="main-tooltip"
          />
        </span>
        <PresetForm>
          <label htmlFor="preset-select">Preset date range:</label>
          <select
            id="preset-select"
            value={DEFAULT_PRESET_OPTION}
            onChange={handlePresetChange}
          >
            {presetDateOptions}
          </select>
        </PresetForm>
        <DateForm>
          <FormColumn>
            <label htmlFor="start">Start date</label>

            <input
              type="date"
              id="start"
              name="date-range-start"
              value={startDate}
              min={config.min_date}
              max={maxDateISO}
              onChange={handleStartChange}
            />
          </FormColumn>
          <FormColumn>
            <label htmlFor="end">End date</label>

            <input
              type="date"
              id="end"
              name="date-range-end"
              value={endDate}
              min={config.min_date}
              max={maxDateISO}
              onChange={handleEndChange}
            />
          </FormColumn>
        </DateForm>
      </Container>
    );
  }
);
DateSelect.propTypes = {
  startDate: PropTypes.string.isRequired,
  endDate: PropTypes.string.isRequired,
  updateDateRange: PropTypes.func.isRequired,
  title: PropTypes.string.isRequired,
  hintText: PropTypes.string.isRequired,
};

export default DateSelect;
