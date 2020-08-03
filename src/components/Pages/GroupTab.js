import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';

import KBD from '../Common/KBD';
import TabIndicator from '../Common/TabIndicator';
import AccordionTitle from '../Common/AccordionTitle';
import AccordionWrapper from '../Common/AccordionWrapper';

import VegaLegend from '../Vega/VegaLegend';
import VegaStackedBars from '../Vega/GroupStackPlot';
import DataTableContainer from '../Table/DataTableContainer';
import AcknowledgementsTable from '../Table/AcknowledgementsTable';

// import { GROUP_KEYS } from '../../constants/config';

const GroupTabContainer = styled.div`
  padding-top: 10px;
`;

const HelpText = styled.div`
  margin-bottom: 5px;

  font-weight: normal;
  font-size: 0.9em;
  color: #666;
  line-height: normal;
  p {
    margin-top: 0px;
    margin-bottom: 3px;
  }
`;

const GroupTab = observer(({ width }) => {
  const { configStore } = useStores();

  return (
    <GroupTabContainer>
      <AccordionWrapper
        title={
          <AccordionTitle>
            <span>Legend</span>
          </AccordionTitle>
        }
        defaultCollapsed={false}
        maxHeight={'500px'}
      >
        <HelpText>
          <p>
            Items in the legend represent <b>{configStore.getGroupLabel()}s</b>.
            Click to select one, or hold <KBD>Shift</KBD> and click to select
            multiple {configStore.getGroupLabel()}s. Selected{' '}
            {configStore.getGroupLabel()}s will be highlighted in the plot and
            table below, as well as in the{' '}
            <TabIndicator>Compare Locations</TabIndicator> tab.
          </p>
        </HelpText>
        <VegaLegend />
      </AccordionWrapper>
      <AccordionWrapper
        title={
          <AccordionTitle>
            <span>Plot</span>
          </AccordionTitle>
        }
        defaultCollapsed={false}
        maxHeight={'1200px'}
      >
        <HelpText>
          <p>
            The plot shows sequences grouped by their respective{' '}
            <b>{configStore.getGroupLabel()}</b> and plotted over time. Click to
            select one, or hold <KBD>Shift</KBD> and click to select multiple{' '}
            {configStore.getGroupLabel()}s. Selected{' '}
            {configStore.getGroupLabel()}s will be highlighted in the legend and
            table below, as well as in the{' '}
            <TabIndicator>Compare Locations</TabIndicator> tab.
          </p>
        </HelpText>
        <VegaStackedBars width={width - 150} />
      </AccordionWrapper>
      {/*configStore.groupKey === GROUP_KEYS.GROUP_LINEAGE && (
        <VegaTree width={width} data={dataStore.caseDataAggGroup} />
      )*/}
      <AccordionWrapper
        title={
          <AccordionTitle>
            <span>Table</span>
          </AccordionTitle>
        }
        defaultCollapsed={false}
        maxHeight={'1200px'}
      >
        <HelpText>
          <p>
            The table shows the counts and associated mutations of each{' '}
            <b>{configStore.getGroupLabel()}</b>. Click a table row to select
            one, or hold <KBD>Shift</KBD> and click to select multiple{' '}
            {configStore.getGroupLabel()}s. Selected{' '}
            {configStore.getGroupLabel()}s will be highlighted in the legend and
            plot, as well as in the{' '}
            <TabIndicator>Compare Locations</TabIndicator> tab.
          </p>
        </HelpText>
        <DataTableContainer />
      </AccordionWrapper>
      <AccordionWrapper
        title="acknowledgements"
        defaultCollapsed={true}
        maxHeight={'1200px'}
      >
        <AcknowledgementsTable />
      </AccordionWrapper>
    </GroupTabContainer>
  );
});
GroupTab.propTypes = {
  width: PropTypes.number,
};
GroupTab.defaultProps = {
  width: 100,
};

export default GroupTab;
