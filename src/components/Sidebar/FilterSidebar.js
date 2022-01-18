import React from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import Header from './Header';
import TabBar from './TabBar';
import StatusBox from './StatusBox';
import Footer from './Footer';

import { FilterSidebarContainer } from './FilterSidebar.styles';

const FilterSidebar = observer(() => {
  const { UIStore } = useStores();

  const onTabChange = (tab) => {
    UIStore.setActiveTab(tab);
  };

  return (
    <FilterSidebarContainer>
      <Header />
      <TabBar activeTab={UIStore.activeTab} onTabChange={onTabChange} />
      <StatusBox />
      <Footer />
    </FilterSidebarContainer>
  );
});

export default FilterSidebar;
