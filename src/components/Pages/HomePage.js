import React from 'react';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import { connect } from '../../stores/connect';
import { onMobileDevice } from '../../utils/device';

import FilterSidebar from '../Sidebar/FilterSidebar';
import DefaultSidebar from '../Sidebar/DefaultSidebar';
import CGLogo from '../../assets/images/cg_logo_v13.png';

import ExampleTab from './ExampleTab';
import GroupTab from './GroupTab';
import LocationTab from './LocationTab';
import AboutTab from './AboutTab';
import MethodologyTab from './MethodologyTab';
import RelatedProjectsTab from './RelatedProjectsTab';
import SequencingEffortsTab from './SequencingEffortsTab';
import KeyListener from '../KeyListener';

import { TABS } from '../../constants/defs.json';

const HomePageDiv = styled.div`
  display: grid;
  grid-template-columns: [col1] 250px [col2] 180px [col3] auto [col4];
  grid-template-rows: [row1] auto [row2];
  width: 100vw;
  position: relative;
  overflow-y: hidden;
`;

const PlotContainer = styled.div`
  grid-column: ${({ showDefaultSidebar }) =>
    showDefaultSidebar ? 'col2 / col4' : 'col3 / col4'};
  grid-row: row1 / row2;
  display: flex;
  flex-direction: column;
  width: 100%;
  height: 100vh;
  box-sizing: border-box;
  position: relative;
  overflow-y: scroll;
  border-left: 1px #eaeaea solid;
`;

const HomePage = observer(({ UIStore }) => {
  const renderTab = () => {
    if (UIStore.activeTab === TABS.TAB_GROUP) {
      return <GroupTab />;
    } else if (UIStore.activeTab === TABS.TAB_LOCATION) {
      return <LocationTab />;
    } else if (UIStore.activeTab === TABS.TAB_EXAMPLE) {
      return <ExampleTab />;
    } else if (UIStore.activeTab === TABS.TAB_ABOUT) {
      return <AboutTab />;
    } else if (UIStore.activeTab === TABS.TAB_METHODOLOGY) {
      return <MethodologyTab />;
    } else if (UIStore.activeTab === TABS.TAB_RELATED) {
      return <RelatedProjectsTab />;
    } else if (UIStore.activeTab === TABS.TAB_GLOBAL_SEQUENCES) {
      return <SequencingEffortsTab />;
    }
  };

  if (onMobileDevice()) {
    return (
      <div
        style={{
          display: 'flex',
          flexDirection: 'column',
          alignItems: 'center',
          marginTop: 10,
        }}
      >
        <img src={CGLogo}></img>
        <p style={{ margin: '20px' }}>
          COVID-19 CG is designed for large screen devices due to the highly
          detailed analyses presented on the browser. Please view this site on a
          laptop or computer for the best user experience.
        </p>
        <p
          style={{
            margin: '20px',
            marginTop: '10px',
            fontWeight: 'normal',
            fontSize: '0.9em',
            lineHeight: 'normal',
          }}
        >
          COVID-19 CG can also be explored on iPad and larger tablets if browser
          settings are switched to &quot;Request Desktop Version&quot; and the
          device is in landscape mode. Exploring COVID-19 CG on handheld mobile
          devices will result in excessively limited browser functionality.
        </p>
      </div>
    );
  }

  const showDefaultSidebar =
    UIStore.activeTab !== TABS.TAB_GROUP &&
    UIStore.activeTab !== TABS.TAB_LOCATION;

  return (
    <>
      <KeyListener />
      <HomePageDiv>
        {showDefaultSidebar ? <DefaultSidebar /> : <FilterSidebar />}
        <PlotContainer showDefaultSidebar={showDefaultSidebar}>
          {renderTab()}
        </PlotContainer>
      </HomePageDiv>
    </>
  );
});

// eslint-disable-next-line react/display-name
export default connect(HomePage);
