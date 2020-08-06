import React, { useState } from 'react';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import { connect } from '../../stores/connect';
// import _ from 'underscore';
import useDimensions from 'react-use-dimensions';

import ReactTooltip from 'react-tooltip';
import QuestionButton from '../Buttons/QuestionButton';
import Header from '../FilterSidebar/Header';
import GroupBySelect from '../FilterSidebar/GroupBySelect';
import CoordinateSelect from '../FilterSidebar/CoordinateSelect';
import MetaFieldSelect from '../FilterSidebar/MetaFieldSelect';
import DropdownContainer from '../FilterSidebar/DropdownContainer';
import SplashScreenModal from '../Modals/SplashScreenModal';
//import initial_entropy_spec from '../vega/barplot_v3.vl.json';
// import SideBar from './Sidebar';
// import VegaTree from './VegaTree';
import TabBar from '../TabBar';
import SidebarAccordionWrapper from '../FilterSidebar/SidebarAccordionWrapper';

import ExampleTab from './ExampleTab';
import GroupTab from './GroupTab';
import LocationTab from './LocationTab';
import AboutTab from './AboutTab';
import MethodologyTab from './MethodologyTab';
import RelatedProjectsTab from './RelatedProjectsTab';
import Footer from '../Footer';
import FilterDataIntoOther from '../FilterSidebar/FilterDataIntoOther';
import KeyListener from '../KeyListener';

import { TABS } from '../../constants/UI';

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
  overflow-y: hidden;

  .filter-sidebar-tooltip {
    background-color: #fff;
    font-weight: normal;
    p {
      margin-top: 2px;
      margin-bottom: 2px;
    }
  }
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

const HomePage = observer(({ UIStore }) => {
  const [ref, { width }] = useDimensions();

  const [modalIsOpen, setIsOpen] = useState(true);
  const openModal = (e) => {
    if (e !== undefined) {
      e.preventDefault();
    }

    setIsOpen(true);
  };
  const afterOpenModal = () => {
    // references are now sync'd and can be accessed.
    // subtitle.style.color = '#f00';
  };
  const closeModal = () => {
    setIsOpen(false);
  };

  const onTabChange = (tab) => {
    UIStore.setActiveTab(tab);
  };

  const renderTab = () => {
    if (UIStore.activeTab === TABS.TAB_GROUP) {
      return <GroupTab width={width} />;
    } else if (UIStore.activeTab === TABS.TAB_LOCATION) {
      return <LocationTab width={width} />;
    } else if (UIStore.activeTab === TABS.TAB_EXAMPLE) {
      return <ExampleTab width={width} />;
    } else if (UIStore.activeTab === TABS.TAB_ABOUT) {
      return <AboutTab />;
    } else if (UIStore.activeTab === TABS.TAB_METHODOLOGY) {
      return <MethodologyTab />;
    } else if (UIStore.activeTab === TABS.TAB_RELATED) {
      return <RelatedProjectsTab />;
    }
  };

  return (
    <>
      <KeyListener />
      <SplashScreenModal
        isOpen={modalIsOpen}
        onAfterOpen={afterOpenModal}
        onRequestClose={closeModal}
      />
      <HomePageDiv>
        {/* <SideBar /> */}
        <FilterSidebar>
          <ReactTooltip
            className="filter-sidebar-tooltip"
            id="tooltip-filter-sidebar"
            type="light"
            effect="solid"
            border={true}
            borderColor="#888"
          />
          <Header />
          <GroupBySelect />
          <SidebarAccordionWrapper
            title="Collapse low frequency data"
            defaultCollapsed={true}
            maxHeight={'250px'}
          >
            <FilterDataIntoOther />
          </SidebarAccordionWrapper>
          <SidebarAccordionWrapper
            title="Genomic coordinates"
            defaultCollapsed={false}
            maxHeight={'420px'}
          >
            <CoordinateSelect />
          </SidebarAccordionWrapper>
          <SidebarAccordionWrapper
            title={
              <div>
                Filter by metadata (advanced)
                <QuestionButton
                  data-tip='<p>By default, no filtering is applied on sequence metadata (Default is select all)</p><p>Metadata is dependent on the data submitter, so many fields may be missing and marked as "Unknown".</p>'
                  data-html={true}
                  data-for="tooltip-filter-sidebar"
                />
              </div>
            }
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
          <TabBar activeTab={UIStore.activeTab} onTabChange={onTabChange} />
          {renderTab()}
          <Footer openModal={openModal} />
        </PlotContainer>
      </HomePageDiv>
    </>
  );
});

// eslint-disable-next-line react/display-name
export default connect(HomePage);
