import React from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { intToISO } from '../../utils/date';
import { MIN_DATE } from '../../constants/defs.json';

import { Container, DateForm, FormColumn } from './DateSelect.styles';

import ReactTooltip from 'react-tooltip';
import QuestionButton from '../Buttons/QuestionButton';

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
