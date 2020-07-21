import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import styled from 'styled-components';
// import _ from 'underscore';
import useDimensions from 'react-use-dimensions';

import CoordinateSelect from './CoordinateSelect';
import GroupBySelect from './GroupBySelect';
import MetaFieldSelect from './MetaFieldSelect';
import DropdownContainer from './DropdownContainer';

//import initial_entropy_spec from '../vega/barplot_v3.vl.json';

import { connect } from '../stores/connect';
import DataTableContainer from './Table/DataTableContainer';
import Header from './Header';
// import SideBar from './Sidebar';
import { asyncStates } from '../stores/uiStore';
import SkeletonElement from './SkeletonElement';
import LoadingSpinner from './LoadingSpinner';
import VegaLegend from './Vega/VegaLegend';
// import VegaTree from './VegaTree';
import AccordionWrapper from './AccordionWrapper';
import SidebarAccordionWrapper from './SidebarAccordionWrapper';
import VegaStackedBars from './Vega/VegaStackedBars';
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

const Footer = styled.div`
  margin-top: auto;
  display: flex;
  background-color: #f8f8f8;

  margin-left: -10px;
  padding: 5px;
  border-top: 1px solid #ccc;

  font-size: 0.85rem;
`;

const HomePage = observer(({ uiStore }) => {
  const [ref, { width }] = useDimensions();

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
          title="Plot"
          defaultCollapsed={false}
          maxHeight={'1200px'}
        >
          <VegaStackedBars width={width - 150} />
        </AccordionWrapper>
      );
    }
  };

  return (
    <>
      <HomePageDiv>
        {/* <SideBar /> */}
        <FilterSidebar>
          <Header />
          <GroupBySelect />
          <SidebarAccordionWrapper
            title="Genomic coordinates"
            defaultCollapsed={false}
            maxHeight={'420px'}
          >
            <CoordinateSelect />
          </SidebarAccordionWrapper>
          <SidebarAccordionWrapper
            title="Filter sequences by"
            defaultCollapsed={true}
            maxHeight={'220px'}
          >
            <MetaFieldSelect />
          </SidebarAccordionWrapper>

          {/*<SidebarAccordionWrapper
            title="Selected locations"
            defaultCollapsed={false}
          >
            <DropdownContainer />
          </SidebarAccordionWrapper>*/}
          <DropdownContainer />
        </FilterSidebar>

        <PlotContainer ref={ref}>
          <AccordionWrapper
            title="Legend"
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
            title="Table"
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
