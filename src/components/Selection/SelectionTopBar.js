import React, { useState } from 'react';
// import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { ASYNC_STATES } from '../../constants/defs.json';

import SkeletonElement from '../Common/SkeletonElement';
import SelectSequencesModal from '../Modals/SelectSequencesModal';
import GroupBySelect from '../Selection/GroupBySelect';
import DownloadDataButton from '../Sidebar/DownloadDataButton';

import {
  Container,
  SelectSequencesButton,
  StatusBox,
  StatusBlock,
} from './SelectionTopBar.styles';

const SelectionTopBar = observer(() => {
  const { UIStore, configStore, dataStore } = useStores();

  const [modalActive, setModalActive] = useState(false);

  const hideModal = () => {
    // Don't close the modal if we're in the middle of a request
    if (
      UIStore.caseDataState !== ASYNC_STATES.SUCCEEDED &&
      UIStore.caseDataState !== ASYNC_STATES.FAILED
    ) {
      return;
    }

    setModalActive(false);
  };
  const showModal = () => {
    setModalActive(true);
  };

  const onChangeGroupKey = (groupKey) => {
    configStore.applyPendingChanges({ groupKey });
  };
  const onChangeDnaOrAa = (dnaOrAa) => {
    configStore.applyPendingChanges({ dnaOrAa });
  };
  const onReferenceChange = (selectedReference) => {
    configStore.applyPendingChanges({ selectedReference });
  };
  const loading = UIStore.caseDataState === ASYNC_STATES.STARTED;

  let statusBox = (
    <StatusBox>
      <StatusBlock>
        <b>{dataStore.numSequencesAfterAllFiltering}</b> sequences fetched in{' '}
        {dataStore.timeToFetch} s.
      </StatusBlock>
    </StatusBox>
  );
  if (loading) {
    statusBox = (
      <div style={{ width: '150px', padding: '5px' }}>
        <SkeletonElement delay={2} height={40}></SkeletonElement>
      </div>
    );
  }

  return (
    <Container>
      <SelectSequencesButton onClick={showModal} disabled={loading}>
        Filter Sequences
      </SelectSequencesButton>
      <SelectSequencesModal isOpen={modalActive} onRequestClose={hideModal} />

      <GroupBySelect
        groupKey={configStore.groupKey}
        dnaOrAa={configStore.dnaOrAa}
        coordinateMode={configStore.coordinateMode}
        selectedGene={configStore.selectedGene}
        selectedProtein={configStore.selectedProtein}
        selectedReference={configStore.selectedReference}
        onGroupKeyChange={onChangeGroupKey}
        onDnaOrAaChange={onChangeDnaOrAa}
        onReferenceChange={onReferenceChange}
        showExtraGroupText={false}
        disabled={loading}
        direction={'row'}
      />

      {statusBox}

      <div className="spacer"></div>

      <DownloadDataButton disabled={loading} direction={'right'} />
    </Container>
  );
});
SelectionTopBar.propTypes = {};
SelectionTopBar.defaultProps = {};

export default SelectionTopBar;
