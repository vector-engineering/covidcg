import React, { useState } from 'react';
// import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { ASYNC_STATES, COORDINATE_MODES } from '../../constants/defs.json';

import SkeletonElement from '../Common/SkeletonElement';
import SelectSequencesModal from '../Modals/SelectSequencesModal';
import GroupBySelect from '../Selection/GroupBySelect';
import DownloadDataButton from '../Sidebar/DownloadDataButton';

import { getGene, getProtein } from '../../utils/gene_protein';

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
    // Update the selected gene/protein object
    const selectedGene = getGene(
      configStore.selectedGene.name,
      selectedReference
    );
    const selectedProtein = getProtein(
      configStore.selectedProtein.name,
      selectedReference
    );

    // Update residue coordinates
    let residueCoordinates;
    if (configStore.coordinateMode === COORDINATE_MODES.COORD_GENE) {
      residueCoordinates = [[1, selectedGene.len_aa]];
    } else if (configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
      residueCoordinates = [[1, selectedProtein.len_aa]];
    }

    configStore.applyPendingChanges({
      selectedReference,
      selectedGene,
      selectedProtein,
      residueCoordinates,
    });
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
        referenceSelectMaxWidth="100px"
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
