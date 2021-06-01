import React from 'react';
import { observer } from 'mobx-react';

import { useStores } from '../../stores/connect';

import Header from './Header';
import TabIndicator from '../Common/TabIndicator';
import TabBar from './TabBar';
import Footer from './Footer';

import { TABS } from '../../constants/defs.json';

import { SidebarContainer, SidebarChunk } from './DefaultSidebar.styles';

const DefaultSidebar = observer(() => {
  const { configStore, UIStore } = useStores();

  const onTabChange = (tab) => {
    UIStore.setActiveTab(tab);
  };

  return (
    <SidebarContainer>
      <Header />
      <TabBar activeTab={UIStore.activeTab} onTabChange={onTabChange} />
      <SidebarChunk>
        To begin analyzing and visualizing data, select the{' '}
        <TabIndicator tab={TABS.TAB_COMPARE_GROUPS}>
          Compare {configStore.getGroupLabel()}s
        </TabIndicator>{' '}
        tab,
        <br />
        or the{' '}
        <TabIndicator tab={TABS.TAB_COMPARE_LOCATIONS}>
          Compare Locations
        </TabIndicator>{' '}
        tab
      </SidebarChunk>
      <Footer />
    </SidebarContainer>
  );
});

export default DefaultSidebar;
