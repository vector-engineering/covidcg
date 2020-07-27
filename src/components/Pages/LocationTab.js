import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';
import { observer } from 'mobx-react';
import { asyncStates } from '../../stores/uiStore';

import ReactTooltip from 'react-tooltip';
import AccordionTitle from '../Common/AccordionTitle';
import AccordionWrapper from '../Common/AccordionWrapper';
import SkeletonElement from '../Common/SkeletonElement';
import LoadingSpinner from '../Common/LoadingSpinner';

import VegaLegend from '../Vega/VegaLegend';
import LocationGroupPlot from '../Vega/LocationGroupPlot';
import LocationDatePlot from '../Vega/LocationDatePlot';

const LocationTabContainer = styled.div``;

const LocationTab = observer(({ width }) => {
  const { covidStore, uiStore } = useStores();

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
              <span>Location-Date Plot</span>
              <span
                className="question-button"
                data-tip="<p>Click on a legend item, line, or dot to select the location.</p><p>Hold the <kbd>Shift</kbd> key and click to select multiple locations.</p>"
                data-html="true"
                data-for="tooltip-location-group"
              >
                ?
              </span>
            </AccordionTitle>
          }
          defaultCollapsed={false}
          maxHeight={'800px'}
        >
          <LocationDatePlot width={width - 200} />
          <ReactTooltip
            id="tooltip-location-date"
            type="light"
            effect="solid"
            border={true}
            borderColor="#888"
          />
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
              <span>Location-Group Plot</span>
              <span
                className="question-button"
                data-tip="<p>Click on a the location names on the y-axis to select locations.</p><p>Hold the <kbd>Shift</kbd> key and click to select multiple locations.</p><p>Click on a bar to select a lineage/clade/SNP.</p><p>Hold the <kbd>Shift</kbd> key and click to select multiple lineages/clades/SNPs.</p>"
                data-html="true"
                data-for="tooltip-location-group"
              >
                ?
              </span>
            </AccordionTitle>
          }
          defaultCollapsed={false}
          maxHeight={'500px'}
        >
          <LocationGroupPlot width={width - 300} />
          <ReactTooltip
            id="tooltip-location-group"
            type="light"
            effect="solid"
            border={true}
            borderColor="#888"
          />
        </AccordionWrapper>
      );
    }
  };

  return (
    <LocationTabContainer>
      <ReactTooltip
        id="tooltip-locations"
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
              data-for="tooltip-locations"
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
