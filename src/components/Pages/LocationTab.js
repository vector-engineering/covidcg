import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';
import { observer } from 'mobx-react';
import { asyncStates } from '../../stores/uiStore';

import AccordionWrapper from '../AccordionWrapper';
import SkeletonElement from '../SkeletonElement';
import LoadingSpinner from '../LoadingSpinner';

import VegaLegend from '../Vega/VegaLegend';
import LocationGroupPlot from '../Vega/LocationGroupPlot';
import LocationDatePlot from '../Vega/LocationDatePlot';

const LocationTabContainer = styled.div``;

const SelectedGroupsContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;

  background-color: #f8f8f8;
  margin-left: 20px;
  border: 1px solid #aaa;
  margin-right: 20px;
  margin-bottom: 10px;
  border-radius: 2px;
  padding: 10px;

  .selected-groups-title {
    margin-right: 10px;
  }

  .group-list {
    display: flex;
    flex-direction: row;
    align-items: center;
    flex-flow: wrap;
    justify-content: start;
  }
`;

const GroupItem = styled.div`
  background-color: #fff;
  border: 1px solid #ddd;
  padding: 3px 8px;
  margin-right: 5px;
`;

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
          <LocationGroupPlot width={width - 300} />
        </AccordionWrapper>
      );
    }
  };

  const groupElements = [];
  covidStore.selectedGroups.forEach((group) => {
    groupElements.push(<GroupItem key={group.group}>{group.group}</GroupItem>);
  });
  let selectedGroupsTitle = 'Selected ';
  if (covidStore.groupKey === 'lineage') {
    selectedGroupsTitle += 'Lineages';
  } else if (covidStore.groupKey === 'clade') {
    selectedGroupsTitle += 'Clades';
  } else if (covidStore.groupKey === 'snp') {
    if (covidStore.dnaOrAa === 'dna') {
      selectedGroupsTitle += 'NT SNPs';
    } else {
      selectedGroupsTitle += 'AA SNPs';
    }
  }

  return (
    <LocationTabContainer>
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
      {/* <SelectedGroupsContainer>
        <span className="selected-groups-title">{selectedGroupsTitle}:</span>
        <div className="group-list">{groupElements}</div>
      </SelectedGroupsContainer> */}

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
