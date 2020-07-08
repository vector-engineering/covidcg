import React, { useState } from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import _ from 'underscore';
import useDimensions from 'react-use-dimensions';

import GeneSelect from './GeneSelect';
import GroupBySelect from './GroupBySelect';
import DropdownContainer from './DropdownContainer';

//import initial_entropy_spec from '../vega/barplot_v3.vl.json';
import areaStackSpecInitial from '../vega/area_stack.vl.json';

import { connect } from '../stores/connect';
import DataTableContainer from './Table/DataTableContainer';
import Header from './Header';
import SideBar from './Sidebar';
import { asyncStates } from '../stores/uiStore';
import SkeletonElement from './SkeletonElement';
import LoadingSpinner from './LoadingSpinner';
import VegaLegend from './VegaLegend';
import VegaTree from './VegaTree';
import AccordionWrapper from './AccordionWrapper';
import VegaWrapper from './VegaWrapper';
import AcknowledgementsTable from './AcknowledgementsTable';

const HomePageDiv = styled.div`
  display: grid;
  grid-template-columns: [col1] 300px [col2] calc(100vw - 300px) [col3];
  grid-template-rows: [row1] auto [row2];
  height: 100vh;
  width: 100vw;
  position: relative;
  overflow-y: hidden;
`;
const FilterSidebar = styled.div`
  grid-column: col1 / col2;
  grid-row: row1 / row2;

  background-color: #f8f8f8;
  //padding-right: 10px;
  padding-bottom: 15px;
  border-right: 1px solid #aaa;
  display: flex;
  flex-direction: column;
  height: 100vh;
`;
const PlotContainer = styled.div`
  grid-column: col2 / col3;
  grid-row: row1 / row2;

  display: flex;
  flex-direction: column;
  width: 100%;
  max-height: 100vh;
  box-sizing: border-box;

  padding-left: 10px;
  padding-top: 10px;

  position: relative;

  overflow-y: scroll;

  .vega-embed {
    width: calc(100% - 110px);
  }
`;
const PlotOptions = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-start;

  .area-stack-title {
    font-size: 1.25em;
    margin-right: 10px;
    padding-right: 10px;
    padding-left: 18px;

    border-right: 1px solid #ccc;
  }
`;

const AreaStackSelectContainer = styled.div`
  font-weight: normal;
  select {
    margin-left: 0.65em;
    padding: 1px 4px;
    border-radius: 3px;
  }
`;

const Footer = styled.div`
  margin-top: auto;
  display: flex;
  background-color: #f8f8f8;

  margin-left: -10px;
  padding: 5px;
  border-top: 1px solid #ccc;

  font-size: 0.85rem;
`;

const AreaStackModeSelect = ({ mode, onChange }) => {
  return (
    <AreaStackSelectContainer>
      <label>
        Display mode:
        <select value={mode} onChange={onChange}>
          <option value="counts">Counts</option>
          <option value="percentages">Percentages</option>
        </select>
      </label>
    </AreaStackSelectContainer>
  );
};
AreaStackModeSelect.propTypes = {
  mode: PropTypes.string.isRequired,
  onChange: PropTypes.func.isRequired,
};

const HomePage = observer(({ covidStore, uiStore }) => {
  // 'percentages' or 'counts'
  const [areaStackMode, setAreaStackMode] = useState('percentages');
  const [ref, { width }] = useDimensions();

  const handleGroupingChange = (groupKey, dnaOrAa) => {
    covidStore.changeGrouping(groupKey, dnaOrAa);
  };

  const handleGeneChange = (gene) => {
    covidStore.selectGene(gene);
  };

  const handleBrush = (...args) => {
    //console.log(args);
    // this.setState({
    //   info: JSON.stringify(args),
    // });
    covidStore.selectDateRange(
      Object.prototype.hasOwnProperty.call(args[1], 'date')
        ? args[1].date
        : [-1, -1]
    );
  };

  const onChangeAreaStackMode = (event) => setAreaStackMode(event.target.value);

  // Make a deep copy of the Vega spec so we can edit it
  const areaStackSpec = JSON.parse(JSON.stringify(areaStackSpecInitial));

  if (areaStackMode === 'percentages') {
    areaStackSpec['vconcat'][0]['encoding']['y']['stack'] = 'normalize';
  }

  // Adapt labels to groupings
  if (covidStore.groupKey === 'lineage') {
    // y-axis title
    areaStackSpec['vconcat'][0]['encoding']['y']['axis']['title'] =
      (areaStackMode === 'percentages' ? 'Percent ' : '') +
      'Sequences by Lineage';
    // Tooltip title
    areaStackSpec['vconcat'][0]['encoding']['tooltip'][0]['title'] = 'Lineage';
  } else if (covidStore.groupKey === 'snp') {
    // y-axis title
    areaStackSpec['vconcat'][0]['encoding']['y']['axis']['title'] =
      (areaStackMode === 'percentages' ? 'Percent ' : '') +
      'Sequences by ' +
      (covidStore.dnaOrAa === 'dna' ? 'NT' : 'AA') +
      ' SNP';
    // Tooltip title
    areaStackSpec['vconcat'][0]['encoding']['tooltip'][0]['title'] =
      (covidStore.dnaOrAa === 'dna' ? 'NT' : 'AA') + ' SNP';
  }

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
          <SkeletonElement delay={2} height={'400px'}>
            <LoadingSpinner />
          </SkeletonElement>
        </div>
      );
    } else {
      return (
        <AccordionWrapper
          title="plot"
          defaultCollapsed={false}
          maxHeight={'1200px'}
        >
          <div style={{ width: `${width}px` }}>
            <VegaWrapper
              data={{
                case_data: covidStore.caseData,
              }}
              spec={areaStackSpec}
              signalListeners={{
                brush: _.debounce(handleBrush, 500),
              }}
            />
          </div>
        </AccordionWrapper>
      );
    }
  };

  let areaStackTitle = 'Lineage ';
  if (covidStore.groupKey === 'lineage') {
    areaStackTitle = 'Lineage ';
  } else if (covidStore.groupKey === 'snp') {
    areaStackTitle = 'SNP ';
  }
  areaStackTitle += areaStackMode === 'percentages' ? 'Percentages' : 'Counts';
  areaStackTitle += ' Over Time';

  return (
    <>
      <HomePageDiv>
        {/* <SideBar /> */}
        <FilterSidebar>
          <Header />
          <GroupBySelect
            groupKey={covidStore.groupKey}
            dnaOrAa={covidStore.dnaOrAa}
            onChange={handleGroupingChange}
          />
          <GeneSelect
            genes={covidStore.genes}
            selectedGene={covidStore.selectedGene}
            onChange={handleGeneChange}
          />

          <DropdownContainer />
        </FilterSidebar>

        <PlotContainer ref={ref}>
          <PlotOptions>
            <span className="area-stack-title">{areaStackTitle}</span>
            <AreaStackModeSelect
              mode={areaStackMode}
              onChange={onChangeAreaStackMode}
            />
          </PlotOptions>
          <br />
          <AccordionWrapper
            title="legend"
            defaultCollapsed={false}
            maxHeight={'500px'}
          >
            <VegaLegend />
          </AccordionWrapper>
          {renderPlotContent()}
          {covidStore.groupKey === 'lineage' && (
            <VegaTree width={width} data={covidStore.caseDataAggGroup} />
          )}

          <AccordionWrapper
            title="table"
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

          <Footer>
            <div className="gisaid-daa">
              Data use subject to the{' '}
              <a
                href="https://www.gisaid.org/"
                target="_blank"
                rel="noopener noreferrer"
              >
                GISAID
              </a>{' '}
              EpiCov™{' '}
              <a
                href="https://www.gisaid.org/registration/terms-of-use/"
                target="_blank"
                rel="noopener noreferrer"
              >
                Database Access Agreement
              </a>
            </div>
          </Footer>
        </PlotContainer>
      </HomePageDiv>
    </>
  );
});

HomePage.propTypes = {
  covidStore: PropTypes.object.isRequired,
  router: PropTypes.object.isRequired,
};

// eslint-disable-next-line react/display-name
export default connect(HomePage);
