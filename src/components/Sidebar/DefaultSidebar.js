import React from 'react';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';

import Header from './Header';
import TabIndicator from '../Common/TabIndicator';

import { TABS } from '../../constants/UI';

const SidebarContainer = styled.div`
  position: fixed;
  top: 0;
  width: 299px;

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

const SidebarChunk = styled.div`
  margin: 5px 12px;
  padding: 5px 0px;
  font-weight: normal;
`;

const DefaultSidebar = observer(() => {
  const { configStore } = useStores();

  return (
    <SidebarContainer>
      <Header />
      <SidebarChunk>
        To begin analyzing and visualizing data, select the{' '}
        <TabIndicator tab={TABS.TAB_GROUP}>
          Compare {configStore.getGroupLabel()}s
        </TabIndicator>{' '}
        tab,
        <br />
        or the{' '}
        <TabIndicator tab={TABS.TAB_LOCATION}>
          Compare Locations
        </TabIndicator>{' '}
        tab
      </SidebarChunk>
    </SidebarContainer>
  );
});

export default DefaultSidebar;
