import React, { useState, useEffect } from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { intToISO, ISOToInt } from '../../utils/date';
import { MIN_DATE } from '../../constants/defs.json';

import { Container, DateForm, FormColumn } from './DateSelect.styles';

import ReactTooltip from 'react-tooltip';
import QuestionButton from '../Buttons/QuestionButton';

const DateSelect = observer(() => {
  const { configStore } = useStores();

  const minDateISO = intToISO(MIN_DATE);
  const maxDate = new Date().getTime();
  const maxDateISO = intToISO(maxDate);

  const [state, setState] = useState({
    startDate: intToISO(configStore.startDate),
    endDate: intToISO(configStore.endDate),
  });

  // Update state from store
  useEffect(() => {
    setState({
      ...state,
      startDate: intToISO(configStore.startDate),
    });
  }, [configStore.startDate]);
  useEffect(() => {
    setState({
      ...state,
      endDate: intToISO(configStore.endDate),
    });
  }, [configStore.endDate]);

  const updateDateRange = ({ startDate, endDate }) => {
    const validDateRange = ISOToInt(startDate) < ISOToInt(endDate);

    configStore.updateValidDateRange(validDateRange);

    if (validDateRange) {
      configStore.updateStartDate(startDate);
    }
  };

  const handleStartChange = (event) => {
    setState({
      ...state,
      startDate: event.target.value,
    });

    updateDateRange({
      startDate: event.target.value,
      endDate: state.endDate,
    });
  };

  const handleEndChange = (event) => {
    setState({
      ...state,
      endDate: event.target.value,
    });

    updateDateRange({
      startDate: state.startDate,
      endDate: event.target.value,
    });
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
            value={state.startDate}
            min={minDateISO}
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
            value={state.endDate}
            min={minDateISO}
            max={maxDateISO}
            onChange={handleEndChange}
          />
        </FormColumn>
      </DateForm>
    </Container>
  );
});

export default DateSelect;
