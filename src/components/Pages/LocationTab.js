import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';
import { observer } from 'mobx-react';

import ExternalLink from '../Common/ExternalLink';
import KBD from '../Common/KBD';
import TabIndicator from '../Common/TabIndicator';
import SelectBoxText from '../Common/SelectBoxText';
import AccordionWrapper from '../Common/AccordionWrapper';

import VegaLegend from '../Vega/VegaLegend';
import LocationGroupPlot from '../Vega/LocationGroupPlot';
import LocationDatePlot from '../Vega/LocationDatePlot';

import { appConfig, GROUP_COLS, GROUP_SNV } from '../../constants/config';
import { TABS } from '../../constants/UI';

const LocationTabContainer = styled.div`
  padding-top: 10px;
`;

const LocationTab = observer(({ width }) => {
  const { configStore } = useStores();

  return (
    <LocationTabContainer>
      <AccordionWrapper
        title="Legend"
        defaultCollapsed={false}
        maxHeight={'500px'}
        helpText={
          <ul>
            <li>
              Items in the legend represent{' '}
              <b>{configStore.getGroupLabel()}s</b>.
            </li>
            <li>
              Click to select one, or hold <KBD>Shift</KBD> and click to select
              multiple {configStore.getGroupLabel()}s.
            </li>
            <li>
              Selected {configStore.getGroupLabel()}s will be highlighted in the
              plots and table below, as well as in the{' '}
              <TabIndicator tab={TABS.TAB_GROUP}>
                Compare {configStore.getGroupLabel()}s
              </TabIndicator>{' '}
              tab.
              {GROUP_COLS.includes(configStore.groupKey) && (
                <>
                  { appConfig.group_defs[configStore.groupKey].description }
                  <ExternalLink href={appConfig.group_defs[configStore.groupKey].link.href}>
                    {appConfig.group_defs[configStore.groupKey].link.title}
                  </ExternalLink>
                </>
              )}
            </li>
          </ul>
        }
      >
        <VegaLegend />
      </AccordionWrapper>
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
            {configStore.groupKey === GROUP_SNV && (
              <li>
                In SNV mode, matching sequences are defined as sequences having
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
            {configStore.groupKey === GROUP_SNV && (
              <li>
                In SNV mode, sequences are shown as counts and not proportions.
                This is because sequences can be counted multiple times, as a
                single sequence may have multiple SNVs. The relevant data is the
                relative difference between the SNVs – the height of all bars is
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
        <LocationGroupPlot width={width - 300} />
      </AccordionWrapper>
    </LocationTabContainer>
  );
});
LocationTab.propTypes = {
  width: PropTypes.number,
};
LocationTab.defaultProps = {
  width: 100,
};

export default LocationTab;
