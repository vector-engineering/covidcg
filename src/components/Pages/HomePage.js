import React, { useState } from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import styled from 'styled-components';
// import _ from 'underscore';
import useDimensions from 'react-use-dimensions';

import CoordinateSelect from '../FilterSidebar/CoordinateSelect';
import GroupBySelect from '../FilterSidebar/GroupBySelect';
import MetaFieldSelect from '../FilterSidebar/MetaFieldSelect';
import DropdownContainer from '../FilterSidebar/DropdownContainer';
import Modal from 'react-modal';

//import initial_entropy_spec from '../vega/barplot_v3.vl.json';

import { connect } from '../../stores/connect';
import DataTableContainer from '../Table/DataTableContainer';
import Header from '../Header';
// import SideBar from './Sidebar';
import { asyncStates } from '../../stores/uiStore';
import SkeletonElement from '../SkeletonElement';
import LoadingSpinner from '../LoadingSpinner';
import VegaLegend from '../Vega/VegaLegend';
// import VegaTree from './VegaTree';
import StatusBar from '../StatusBar';
import AccordionWrapper from '../AccordionWrapper';
import ReactTooltip from 'react-tooltip';
import SidebarAccordionWrapper from '../LiteMol/SidebarAccordionWrapper';
import VegaStackedBars from '../Vega/VegaStackedBars';
import AcknowledgementsTable from '../AcknowledgementsTable';

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

  position: relative;

  overflow-y: scroll;
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

Modal.setAppElement('#app');

const HomePage = observer(({ uiStore }) => {
  const [ref, { width }] = useDimensions();

  const [modalIsOpen, setIsOpen] = useState(false);
  const openModal = () => {
    setIsOpen(true);
  };
  const afterOpenModal = () => {
    // references are now sync'd and can be accessed.
    // subtitle.style.color = '#f00';
  };
  const closeModal = () => {
    setIsOpen(false);
  };

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
    <>
      <Modal
        isOpen={modalIsOpen}
        onAfterOpen={afterOpenModal}
        onRequestClose={closeModal}
        style={{
          content: {
            top: '50%',
            left: '50%',
            right: 'auto',
            bottom: 'auto',
            marginRight: '-50%',
            transform: 'translate(-50%, -50%)',
          },
        }}
        contentLabel="Example Modal"
      >
        <h2>Hello</h2>
        <button onClick={closeModal}>close</button>
        <div>I am a modal</div>
        <form>
          <input />
          <button>tab navigation</button>
          <button>stays</button>
          <button>inside</button>
          <button>the modal</button>
        </form>
      </Modal>
      <HomePageDiv>
        <ReactTooltip
          id="tooltip-home"
          type="light"
          effect="solid"
          border={true}
          borderColor="#888"
        />
        {/* <SideBar /> */}
        <FilterSidebar>
          <button onClick={openModal}>Open Modal</button>
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
            maxHeight={'240px'}
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
          <StatusBar />
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
