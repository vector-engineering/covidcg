import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';
import { asyncStates } from '../../stores/uiStore';

import ReactTooltip from 'react-tooltip';
import AccordionTitle from '../Common/AccordionTitle';
import AccordionWrapper from '../Common/AccordionWrapper';
import SkeletonElement from '../Common/SkeletonElement';
import LoadingSpinner from '../Common/LoadingSpinner';

import VegaLegend from '../Vega/VegaLegend';
import VegaStackedBars from '../Vega/VegaStackedBars';
import DataTableContainer from '../Table/DataTableContainer';
import AcknowledgementsTable from '../Table/AcknowledgementsTable';

const GroupTabContainer = styled.div``;

const GroupTab = observer(({ width }) => {
  const { uiStore } = useStores();

  const renderPlotContent = () => {
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
              <span
                className="question-button"
                data-tip="<p>Hover over a bar in the plot to highlight it in the legend and table.</p><p>Click on a bar to select it, here and in the legend and table.</p><p>Hold the <kbd>Shift</kbd> key and click to select multiple bars.</p>"
                data-html="true"
                data-for="tooltip-plot"
              >
                ?
              </span>
            </AccordionTitle>
          }
          defaultCollapsed={false}
          maxHeight={'1200px'}
        >
          <ReactTooltip
            id="tooltip-plot"
            type="light"
            effect="solid"
            border={true}
            borderColor="#888"
          />
          <VegaStackedBars width={width - 150} />
        </AccordionWrapper>
      );
    }
  };

  return (
    <GroupTabContainer>
      <ReactTooltip
        id="tooltip-home"
        type="light"
        effect="solid"
        border={true}
        borderColor="#888"
      />
      <AccordionWrapper
        title={
          <AccordionTitle>
            <span>Legend</span>
            <span
              className="question-button"
              data-tip="<p>Hover over an item in the legend to highlight it in the plot and table.</p><p>Click on a legend item to select it, here and in the plot and table.</p><p>Hold the <kbd>Shift</kbd> key and click to select multiple items.</p>"
              data-html="true"
              data-for="tooltip-home"
            >
              ?
            </span>
          </AccordionTitle>
        }
        defaultCollapsed={false}
        maxHeight={'500px'}
      >
        <VegaLegend />
      </AccordionWrapper>
      {renderPlotContent()}
      {/*covidStore.groupKey === 'lineage' && (
        <VegaTree width={width} data={covidStore.caseDataAggGroup} />
      )*/}
      <AccordionWrapper
        title={
          <AccordionTitle>
            <span>Table</span>
            <span
              className="question-button"
              data-tip="<p>Hover over an row in the table to highlight it in the plot and legend.</p><p>Click on a row to select it, here and in the plot and legend.</p><p>Hold the <kbd>Shift</kbd> key and click to select multiple rows.</p>"
              data-html="true"
              data-for="tooltip-home"
            >
              ?
            </span>
          </AccordionTitle>
        }
        defaultCollapsed={false}
        maxHeight={'1200px'}
      >
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
