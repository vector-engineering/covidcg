import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';
import { observer } from 'mobx-react';
import { asyncStates } from '../../stores/uiStore';

import AccordionWrapper from '../AccordionWrapper';
import SkeletonElement from '../SkeletonElement';
import LoadingSpinner from '../LoadingSpinner';

import LocationGroupPlot from '../Vega/LocationGroupPlot';
import LocationDatePlot from '../Vega/LocationDatePlot';

const LocationTabContainer = styled.div``;

const AccordionTitle = styled.span`
  .question-button {
    font-family: monospace;
    font-size: 1em;
    line-height: normal;

    margin-left: 8px;
    padding: 2px 5px;
    background-color: #fff;
    border: 1px solid #ccc;
    border-radius: 2px;
    &:hover {
      background-color: #f8f8f8;
    }
  }
`;

const LocationTab = observer(({ width }) => {
  const { uiStore } = useStores();

  const renderLocationDatePlot = () => {
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
          <SkeletonElement delay={2} height={300}>
            <LoadingSpinner />
          </SkeletonElement>
        </div>
      );
    } else {
      return (
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
          <LocationDatePlot width={width - 200} />
        </AccordionWrapper>
      );
    }
  };

  const renderLocationGroupPlot = () => {
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
          <SkeletonElement delay={2} height={300}>
            <LoadingSpinner />
          </SkeletonElement>
        </div>
      );
    } else {
      return (
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
          <LocationGroupPlot width={width - 200} />
        </AccordionWrapper>
      );
    }
  };

  return (
    <LocationTabContainer>
      {renderLocationDatePlot()}
      {renderLocationGroupPlot()}
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
