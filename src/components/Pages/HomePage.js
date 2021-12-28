import React, { useState, useEffect } from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { onMobileDevice } from '../../utils/device';

import ReactTooltip from 'react-tooltip';
import FilterSidebar from '../Sidebar/FilterSidebar';
import DefaultSidebar from '../Sidebar/DefaultSidebar';
import CGLogo from '../../assets/images/cg_logo_v13.png';

import { TABS, ASYNC_STATES } from '../../constants/defs.json';
import KeyListener from '../KeyListener';
import AsyncErrorModal from '../Modals/AsyncErrorModal';

const CompareGroupsTab = React.lazy(() => import('./CompareGroupsTab'));
const HomeTab = React.lazy(() => import('./HomeTab'));
const CompareLocationsTab = React.lazy(() => import('./CompareLocationsTab'));
const GroupReportTab = React.lazy(() => import('./GroupReportTab'));
const AboutTab = React.lazy(() => import('./AboutTab'));
const MethodologyTab = React.lazy(() => import('./MethodologyTab'));
const RelatedProjectsTab = React.lazy(() => import('./RelatedProjectsTab'));
const SequencingEffortsTab = React.lazy(() => import('./SequencingEffortsTab'));

import { HomePageDiv, PlotContainer } from './HomePage.styles';
import MobileHomePage from '../Mobile/MobileHomePage';

const HomePage = observer(() => {
  const { UIStore } = useStores();
  const [showAsyncError, setShowAsyncError] = useState(false);

  const showFetchErrorModal = () => {
    setShowAsyncError(true);
  };
  const hideFetchErrorModal = () => {
    setShowAsyncError(false);
  };

  useEffect(() => {
    if (UIStore.caseDataState === ASYNC_STATES.FAILED) {
      showFetchErrorModal();
    }
  }, [UIStore.caseDataState]);

  const renderTab = () => {
    if (UIStore.activeTab === TABS.TAB_COMPARE_GROUPS) {
      return <CompareGroupsTab />;
    } else if (UIStore.activeTab === TABS.TAB_COMPARE_LOCATIONS) {
      return <CompareLocationsTab />;
    } else if (UIStore.activeTab === TABS.TAB_GROUP_REPORT) {
      return <GroupReportTab />;
    } else if (UIStore.activeTab === TABS.TAB_EXAMPLE) {
      return <HomeTab />;
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
      <MobileHomePage></MobileHomePage>
      // <div
      //   style={{
      //     display: 'flex',
      //     flexDirection: 'column',
      //     alignItems: 'center',
      //     marginTop: 10,
      //   }}
      // >
      //   <img src={CGLogo}></img>
      //   <p style={{ margin: '20px' }}>
      //     COVID-19 CG is designed for large screen devices due to the highly
      //     detailed analyses presented on the browser. Please view this site on a
      //     laptop or computer for the best user experience.
      //   </p>
      //   <p
      //     style={{
      //       margin: '20px',
      //       marginTop: '10px',
      //       fontWeight: 'normal',
      //       fontSize: '0.9em',
      //       lineHeight: 'normal',
      //     }}
      //   >
      //     COVID-19 CG can also be explored on iPad and larger tablets if browser
      //     settings are switched to &quot;Request Desktop Version&quot; and the
      //     device is in landscape mode. Exploring COVID-19 CG on handheld mobile
      //     devices will result in excessively limited browser functionality.
      //   </p>
      // </div>
    );
  }

  const showDefaultSidebar =
    UIStore.activeTab !== TABS.TAB_COMPARE_GROUPS &&
    UIStore.activeTab !== TABS.TAB_COMPARE_LOCATIONS;

  return (
    <>
      <KeyListener />
      <AsyncErrorModal
        isOpen={showAsyncError}
        onAfterOpen={() => {}}
        onRequestClose={hideFetchErrorModal}
      />
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
          <React.Suspense fallback={<div />}>
            {renderTab()}
          </React.Suspense>
        </PlotContainer>
      </HomePageDiv>
    </>
  );
});

// eslint-disable-next-line react/display-name
export default HomePage;
