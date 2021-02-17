import React from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { intToISO } from '../../utils/date';
import { MIN_DATE } from '../../constants/defs.json';

import {
  Container,
  DateForm,
  FormColumn,
  PresetForm,
} from './DateSelect.styles';

import ReactTooltip from 'react-tooltip';
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
    Since pandemic start
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

const DateSelect = observer(() => {
  const { configStore } = useStores();

  const maxDate = new Date().getTime();
  const maxDateISO = intToISO(maxDate);

  const handleStartChange = (event) => {
    configStore.updateStartDate(event.target.value);
  };

  const handleEndChange = (event) => {
    configStore.updateEndDate(event.target.value);
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
      startDate = MIN_DATE;
    }

    configStore.updateStartDate(startDate);
    configStore.updateEndDate(intToISO(endDate));
  };

  return (
    <Container>
      <ReactTooltip
        id="date-select-tooltip"
        type="light"
        effect="solid"
        border={true}
        borderColor="#888"
      />
      <span className="title">
        Date Range
        <QuestionButton
          data-tip={`<p>Filter for sequences that were collected between the two dates.<br/>Dates are start- and end-inclusive, i.e., [start, end]</p>`}
          data-html={true}
          data-for="date-select-tooltip"
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
            value={configStore.startDate}
            min={MIN_DATE}
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
            value={configStore.endDate}
            min={MIN_DATE}
            max={maxDateISO}
            onChange={handleEndChange}
          />
        </FormColumn>
      </DateForm>
    </Container>
  );
});

export default DateSelect;
