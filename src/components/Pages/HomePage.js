import React, { useState } from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import styled from 'styled-components';
// import _ from 'underscore';
import useDimensions from 'react-use-dimensions';
import { connect } from '../../stores/connect';

import Header from '../FilterSidebar/Header';
import GroupBySelect from '../FilterSidebar/GroupBySelect';
import CoordinateSelect from '../FilterSidebar/CoordinateSelect';
import MetaFieldSelect from '../FilterSidebar/MetaFieldSelect';
import DropdownContainer from '../FilterSidebar/DropdownContainer';
import SplashScreenModal from '../Modals/SplashScreenModal';
//import initial_entropy_spec from '../vega/barplot_v3.vl.json';
// import SideBar from './Sidebar';
// import VegaTree from './VegaTree';
import StatusBar from '../StatusBar';
import SidebarAccordionWrapper from '../LiteMol/SidebarAccordionWrapper';

import GroupTab from './GroupTab';
import AboutTab from './AboutTab';
import Footer from '../Footer';

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

const HomePage = observer(() => {
  const [ref, { width }] = useDimensions();

  const [modalIsOpen, setIsOpen] = useState(true);
  const [activeTab, setActiveTab] = useState('group');
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
    setActiveTab(tab);
  };

  const renderTab = () => {
    if (activeTab === 'group') {
      return <GroupTab width={width} />;
    } else if (activeTab === 'location') {
      return <div></div>;
    } else if (activeTab === 'about') {
      return <AboutTab />;
    }
  };

  return (
    <>
      <SplashScreenModal
        isOpen={modalIsOpen}
        onAfterOpen={afterOpenModal}
        onRequestClose={closeModal}
      />
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
          <StatusBar activeTab={activeTab} onTabChange={onTabChange} />
          {renderTab()}
          <Footer openModal={openModal} />
        </PlotContainer>
      </HomePageDiv>
    </>
  );
});

HomePage.propTypes = {
  router: PropTypes.object.isRequired,
};

// eslint-disable-next-line react/display-name
export default connect(HomePage);
