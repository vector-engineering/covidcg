import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';
import { asyncStates } from '../../stores/uiStore';

import KBD from '../Common/KBD';
import TabIndicator from '../Common/TabIndicator';
import AccordionTitle from '../Common/AccordionTitle';
import AccordionWrapper from '../Common/AccordionWrapper';
import SkeletonElement from '../Common/SkeletonElement';
import LoadingSpinner from '../Common/LoadingSpinner';

import VegaLegend from '../Vega/VegaLegend';
import VegaStackedBars from '../Vega/VegaStackedBars';
import DataTableContainer from '../Table/DataTableContainer';
import AcknowledgementsTable from '../Table/AcknowledgementsTable';

const GroupTabContainer = styled.div``;

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

const PlotContent = observer(({ width }) => {
  const { covidStore, uiStore } = useStores();

  if (uiStore.caseDataState === asyncStates.STARTED) {
    return (
      <div
        style={{
          paddingTop: '12px',
          paddingRight: '24px',
          paddingLeft: '12px',
          paddingBottom: '24px',
        }}
      >
        <SkeletonElement delay={2} height={400}>
          <LoadingSpinner />
        </SkeletonElement>
      </div>
    );
  } else {
    return (
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
            <b>{covidStore.getGroupLabel()}</b> and plotted over time. Click to
            select one, or hold <KBD>Shift</KBD> and click to select multiple{' '}
            {covidStore.getGroupLabel()}s. Selected {covidStore.getGroupLabel()}
            s will be highlighted in the legend and table below, as well as in
            the <TabIndicator>Compare Locations</TabIndicator> tab.
          </p>
        </HelpText>
        <VegaStackedBars width={width - 150} />
      </AccordionWrapper>
    );
  }
});

const GroupTab = observer(({ width }) => {
  const { covidStore } = useStores();

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
            Items in the legend represent <b>{covidStore.getGroupLabel()}s</b>.
            Click to select one, or hold <KBD>Shift</KBD> and click to select
            multiple {covidStore.getGroupLabel()}s. Selected{' '}
            {covidStore.getGroupLabel()}s will be highlighted in the plot and
            table below, as well as in the{' '}
            <TabIndicator>Compare Locations</TabIndicator> tab.
          </p>
        </HelpText>
        <VegaLegend />
      </AccordionWrapper>
      <PlotContent width={width} />
      {/*covidStore.groupKey === 'lineage' && (
        <VegaTree width={width} data={covidStore.caseDataAggGroup} />
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
            <b>{covidStore.getGroupLabel()}</b>. Click a table row to select
            one, or hold <KBD>Shift</KBD> and click to select multiple{' '}
            {covidStore.getGroupLabel()}s. Selected {covidStore.getGroupLabel()}
            s will be highlighted in the legend and plot, as well as in the{' '}
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
