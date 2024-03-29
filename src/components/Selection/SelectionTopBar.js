import React, { useState } from 'react';
// import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import {
  ASYNC_STATES,
  COORDINATE_MODES,
  GROUP_MUTATION,
} from '../../constants/defs.json';

import SkeletonElement from '../Common/SkeletonElement';
import SelectSequencesModal from '../Modals/SelectSequencesModal';
import GroupBySelect from '../Selection/GroupBySelect';
import DownloadDataButton from '../Sidebar/DownloadDataButton';

import { getGene, getProtein } from '../../utils/gene_protein';
import { getDefaultReferenceForSubtype } from '../../utils/reference';
import { configStore as initialConfigStore } from '../../constants/initialValues';

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
    // If we just changed from mutation to another grouping,
    // then clear selected group fields
    let selectedGroupFields = configStore.selectedGroupFields;
    if (groupKey !== GROUP_MUTATION && groupKey !== configStore.groupKey) {
      selectedGroupFields = {};
    }
    // If we're switching to mutation mode, go back to
    // default selected group fields
    else if (groupKey === GROUP_MUTATION && groupKey !== configStore.groupKey) {
      selectedGroupFields = initialConfigStore.selectedGroupFields;
    }

    configStore.applyPendingChanges({ groupKey, selectedGroupFields });
  };
  const onChangeDnaOrAa = (dnaOrAa) => {
    configStore.applyPendingChanges({ dnaOrAa });
  };

  const getFeaturesFromNewReference = (newReference) => {
    // Update the selected gene/protein object
    const selectedGene = getGene(configStore.selectedGene.name, newReference);
    const selectedProtein = getProtein(
      configStore.selectedProtein.name,
      newReference
    );

    // Update residue coordinates
    let residueCoordinates;
    if (configStore.coordinateMode === COORDINATE_MODES.COORD_GENE) {
      residueCoordinates = [selectedGene.residue_offset_range];
    } else if (configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
      residueCoordinates = [selectedProtein.residue_offset_range];
    }

    return {
      selectedGene,
      selectedProtein,
      residueCoordinates,
    };
  };

  const onReferenceChange = (selectedReference) => {
    configStore.applyPendingChanges({
      selectedReference,
      ...getFeaturesFromNewReference(selectedReference),
    });
  };
  const onSelectedGroupFieldsChange = (selectedGroupFields) => {
    // Get the first reference from the new subtype
    const selectedReference = getDefaultReferenceForSubtype(
      selectedGroupFields.subtype[0]
    );

    configStore.applyPendingChanges({
      selectedReference,
      selectedGroupFields,
      ...getFeaturesFromNewReference(selectedReference),
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
        selectedGroupFields={configStore.selectedGroupFields}
        onGroupKeyChange={onChangeGroupKey}
        onDnaOrAaChange={onChangeDnaOrAa}
        onReferenceChange={onReferenceChange}
        onSelectedGroupFieldsChange={onSelectedGroupFieldsChange}
        showExtraGroupText={false}
        disabled={loading}
        direction={'row'}
        referenceSelectMaxWidth="120px"
        style={{ flexShrink: 0 }}
      />

      {statusBox}

      <div className="spacer"></div>

      <DownloadDataButton
        disabled={loading}
        direction={'right'}
        style={{ flexShrink: 0 }}
      />
    </Container>
  );
});
SelectionTopBar.propTypes = {};
SelectionTopBar.defaultProps = {};

export default SelectionTopBar;
