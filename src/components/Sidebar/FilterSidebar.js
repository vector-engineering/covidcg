import React, { useState } from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import { ASYNC_STATES } from '../../constants/defs.json';

import Header from './Header';
import Legend from '../Legend/Legend';
import TabBar from './TabBar';
import StatusBox from './StatusBox';
import Footer from './Footer';
import SelectSequencesModal from '../Modals/SelectSequencesModal';

import {
  Container,
  FilterSidebarContainer,
  SelectSequencesButton,
  LegendSidebarContainer,
} from './FilterSidebar.styles';

const FilterSidebar = observer(() => {
  const { UIStore } = useStores();

  const [modalActive, setModalActive] = useState(false);

  const hideModal = () => {
    // Don't close the modal if we're in the middle of a request
    if (UIStore.caseDataState !== ASYNC_STATES.SUCCEEDED) {
      return;
    }

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
        <SelectSequencesButton onClick={showModal}>
          Select Sequences
        </SelectSequencesButton>
        <SelectSequencesModal isOpen={modalActive} onRequestClose={hideModal} />
        <StatusBox />
        <Footer />
      </FilterSidebarContainer>
      <LegendSidebarContainer>
        <Legend />
      </LegendSidebarContainer>
    </Container>
  );
});

export default FilterSidebar;
