import React, { useState } from 'react';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';

import Header from './Header';
import LegendSidebar from './LegendSidebar';
import TabBar from './TabBar';

import Button from '../Buttons/Button';
import SelectSequencesModal from '../Modals/SelectSequencesModal';

const Container = styled.div`
  width: 350px;
  display: flex;
  height: 100vh;
  overflow-y: hidden;
`;

const FilterSidebarContainer = styled.div`
  width: 199px;

  background-color: #f8f8f8;
  //padding-right: 10px;
  padding-bottom: 15px;
  border-right: 1px solid #aaa;
  display: flex;
  flex-direction: column;
  overflow-y: hidden;
  height: 100%;

  .filter-sidebar-tooltip {
    background-color: #fff;
    font-weight: normal;
    p {
      margin-top: 2px;
      margin-bottom: 2px;
    }
  }
`;

const LegendSidebarContainer = styled.div`
  width: 150px;
  height: 100%;
  border-right: 1px solid #aaa;
  padding-bottom: 15px;
`;

const FilterSidebar = observer(() => {
  const { UIStore } = useStores();

  const [modalActive, setModalActive] = useState(false);

  const hideModal = () => {
    setModalActive(false);
  };
  const showModal = () => {
    setModalActive(true);
  };

  const onTabChange = (tab) => {
    UIStore.setActiveTab(tab);
  };

  return (
    <Container>
      <FilterSidebarContainer>
        <Header />
        <TabBar activeTab={UIStore.activeTab} onTabChange={onTabChange} />
        <Button onClick={showModal}>Select Sequences</Button>
        <SelectSequencesModal isOpen={modalActive} onRequestClose={hideModal} />
      </FilterSidebarContainer>
      <LegendSidebarContainer>
        <LegendSidebar />
      </LegendSidebarContainer>
    </Container>
  );
});

export default FilterSidebar;
