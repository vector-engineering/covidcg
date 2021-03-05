import React from 'react';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import { connect } from '../../stores/connect';
import { onMobileDevice } from '../../utils/device';

import ReactTooltip from 'react-tooltip';
import FilterSidebar from '../Sidebar/FilterSidebar';
import DefaultSidebar from '../Sidebar/DefaultSidebar';
import CGLogo from '../../assets/images/cg_logo_v13.png';

import { TABS } from '../../constants/defs.json';
import KeyListener from '../KeyListener';

const GroupTab = React.lazy(() => import('./GroupTab'));
const ExampleTab = React.lazy(() => import('./ExampleTab'));
const LocationTab = React.lazy(() => import('./LocationTab'));
const AboutTab = React.lazy(() => import('./AboutTab'));
const MethodologyTab = React.lazy(() => import('./MethodologyTab'));
const RelatedProjectsTab = React.lazy(() => import('./RelatedProjectsTab'));
const SequencingEffortsTab = React.lazy(() => import('./SequencingEffortsTab'));

const HomePageDiv = styled.div`
  display: grid;
  grid-template-columns: [col1] 250px [col2] 180px [col3] auto [col4];
  grid-template-rows: [row1] auto [row2];
  width: 100vw;
  position: relative;
  overflow-y: hidden;

  .main-tooltip {
    max-width: 300px;

    background-color: #fff;
    font-weight: normal;
    p {
      margin-top: 2px;
      margin-bottom: 2px;
    }
  }
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
        <ReactTooltip
          className="main-tooltip"
          id="main-tooltip"
          type="light"
          effect="solid"
          border={true}
          borderColor="#888"
        />
        {showDefaultSidebar ? <DefaultSidebar /> : <FilterSidebar />}
        <PlotContainer showDefaultSidebar={showDefaultSidebar}>
          <React.Suspense fallback={<div />}>{renderTab()}</React.Suspense>
        </PlotContainer>
      </HomePageDiv>
    </>
  );
});

// eslint-disable-next-line react/display-name
export default connect(HomePage);
