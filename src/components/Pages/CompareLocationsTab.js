import React from 'react';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';
import { observer } from 'mobx-react';
import useDimensions from 'react-use-dimensions';

import KBD from '../Common/KBD';
import SelectBoxText from '../Common/SelectBoxText';
import AccordionWrapper from '../Common/AccordionWrapper';

import LocationGroupPlot from '../Vega/LocationGroupPlot';
import LocationDatePlot from '../Vega/LocationDatePlot';
import NumSeqPerLocationBar from '../Vega/NumSeqPerLocationBar';
import NumSeqPerLocationLine from '../Vega/NumSeqPerLocationLine';

import { GROUP_MUTATION } from '../../constants/defs.json';

const CompareLocationsTabContainer = styled.div`
  padding-top: 10px;
`;

const CompareLocationsTab = observer(() => {
  const { configStore } = useStores();
  const [ref, { width }] = useDimensions();

  const renderLocationDatePlot = () => {
    return (
      <AccordionWrapper
        title="Location-Date Plot"
        defaultCollapsed={false}
        maxHeight={'800px'}
        helpText={
          <ul>
            <li>
              This plot compares the sequence counts or percentages, of the
              selected <b>{configStore.getGroupLabel()}s</b>, between the
              selected locations.
            </li>
            {configStore.groupKey === GROUP_MUTATION && (
              <li>
                In mutation mode, matching sequences are defined as sequences having
                the selected {configStore.getGroupLabel()}s or combination of{' '}
                {configStore.getGroupLabel()}s.
              </li>
            )}
            <li>
              Sequences can be shown as new occurrences{' '}
              <SelectBoxText>New</SelectBoxText>, or as a cumulative count over
              time <SelectBoxText>Cumulative</SelectBoxText>.
            </li>
            <li>
              Sequences can be plotted as raw counts{' '}
              <SelectBoxText>Counts</SelectBoxText>, or as percentages
              normalized to the total number of sequences at that time
              <SelectBoxText>Percentages</SelectBoxText>
            </li>
            <li>
              Sequences are plotted over time by the sample collection date.
              Sequences can be grouped into time bins by{' '}
              <SelectBoxText>Day</SelectBoxText>,{' '}
              <SelectBoxText>Week</SelectBoxText>, or{' '}
              <SelectBoxText>Month</SelectBoxText>
            </li>
            <li>
              <i>Click</i> to highlight one, or hold <KBD>Shift</KBD> and
              <i>click</i> to highlight multiple locations. Highlighted
              locations will be shown in the plot below as well.
            </li>
          </ul>
        }
      >
        <LocationDatePlot width={width - 200} />
      </AccordionWrapper>
    );
  };

  const renderNumSeqPerLocationBarPlot = () => {
    return (
      <AccordionWrapper
        title="Number of Sequences Per Location"
        defaultCollapsed={false}
        maxHeight={'500px'}
        helpText={
          <ul>
            <li>
              This plot shows the number of sequences {' '}
              <b>{configStore.getGroupLabel()}s</b> per location.
            </li>
            <li>
              This can be used to determine biases where large amounts of sequences
              can potentially effect how certain visualizations and data will appear.
            </li>
          </ul>
        }
      >
        <NumSeqPerLocationBar width={width - 200} />
      </AccordionWrapper>
    );
  };

  const renderNumSeqPerLocationLinePlot = () => {
    return (
      <AccordionWrapper
        title="Number of Sequences Per Location Over Time"
        defaultCollapsed={false}
        maxHeight={'500px'}
        helpText={
          <ul>
            <li>
              This plot shows the number of sequences {' '}
              <b>{configStore.getGroupLabel()}s</b> per location over the selected time frames.
            </li>
            <li>
              This can be used to determine biases where large amounts of sequences
              can potentially effect how certain visualizations and data will appear.
            </li>
          </ul>
        }
      >
        <NumSeqPerLocationLine width={width - 200} />
      </AccordionWrapper>
    );
  };

  const renderLocationGroupPlot = () => {
    return (
      <AccordionWrapper
        title="Location-Group Plot"
        defaultCollapsed={false}
        maxHeight={'500px'}
        helpText={
          <ul>
            <li>
              This plot shows the cumulative proportion of{' '}
              <b>{configStore.getGroupLabel()}s</b> per location.
            </li>
            {configStore.groupKey === GROUP_MUTATION && (
              <li>
                In mutation mode, sequences are shown as counts and not proportions.
                This is because sequences can be counted multiple times, as a
                single sequence may have multiple mutations. The relevant data is the
                relative difference between the mutations – the height of all bars is
                not meaningful.
              </li>
            )}
            <li>
              <i>Click</i> to select one, or hold <KBD>Shift</KBD> and{' '}
              <i>click</i> to select multiple {configStore.getGroupLabel()}s.
              Sequences from the selected {configStore.getGroupLabel()}s will be
              shown in the plot above.
            </li>
          </ul>
        }
      >
        <LocationGroupPlot width={width - 200} />
      </AccordionWrapper>
    );
  };

  return (
    <CompareLocationsTabContainer ref={ref}>
      {renderLocationDatePlot()}
      {renderLocationGroupPlot()}
      {renderNumSeqPerLocationBarPlot()}
      {renderNumSeqPerLocationLinePlot()}
    </CompareLocationsTabContainer>
  );
});

export default CompareLocationsTab;
