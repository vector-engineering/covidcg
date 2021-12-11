import React, { useEffect, useRef, useState } from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import PropTypes from 'prop-types';

import { getGene, getProtein } from '../../utils/gene_protein';
import { ISOToInt, intToISO } from '../../utils/date';
import {
  COORDINATE_MODES,
  ASYNC_STATES,
  GROUP_MUTATION,
  DNA_OR_AA,
  MIN_DATE,
} from '../../constants/defs.json';

import Modal from 'react-modal';

import LocationSelect from '../Selection/LocationSelect';
import GroupBySelect from '../Selection/GroupBySelect';
import CoordinateSelect from '../Selection/CoordinateSelect';
import DateSelect from '../Selection/DateSelect';
import MetaFieldSelect from '../Selection/MetaFieldSelect';
import LoadingSpinner from '../Common/LoadingSpinner';

import {
  Overlay,
  ProgressContainer,
  ProgressText,
  TitleContainer,
  HeaderContainer,
  HeaderRow,
  HeaderButtons,
  CancelButton,
  InvalidText,
} from './Modal.styles';

import {
  Wrapper,
  Content,
  Column,
  ApplyButton,
} from './SelectSequencesModal.styles';

Modal.setAppElement('#app');
const NOOP = () => {};

const SelectSequencesContent = observer(({ onRequestClose }) => {
  const { configStore, UIStore, locationDataStore, metadataStore } =
    useStores();
  const sentRequest = useRef(false);

  const [pending, setPending] = useState({
    groupKey: configStore.groupKey,
    dnaOrAa: configStore.dnaOrAa,
    selectedGene: configStore.selectedGene,
    selectedProtein: configStore.selectedProtein,
    selectedPrimers: configStore.selectedPrimers,
    customCoordinates: configStore.customCoordinates,
    validCustomCoordinates: true,
    customSequences: configStore.customSequences,
    validCustomSequences: true,
    residueCoordinates: configStore.residueCoordinates,
    validResidueCoordinates: true,
    coordinateMode: configStore.coordinateMode,
    selectedLocationNodes: configStore.selectedLocationNodes,
    startDate: configStore.startDate,
    endDate: configStore.endDate,
    submStartDate: configStore.submStartDate,
    submEndDate: configStore.submEndDate,
    selectedMetadataFields: configStore.selectedMetadataFields,
    ageRange: configStore.ageRange,
  });

  const changeGrouping = (groupKey, dnaOrAa) => {
    let selectedGene = pending.selectedGene;
    let selectedProtein = pending.selectedProtein;
    // If we switched to non-mutation grouping in AA-mode,
    // then make sure we don't have "All Genes" or "All Proteins" selected
    if (groupKey !== GROUP_MUTATION && dnaOrAa === DNA_OR_AA.AA) {
      if (selectedGene.name === 'All Genes') {
        // Switch back to S gene
        selectedGene = getGene('S');
      }
      if (selectedProtein.name === 'All Proteins') {
        // Switch back to nsp12 protein
        selectedProtein = getProtein('nsp12 - RdRp');
      }
    }

    setPending({
      ...pending,
      groupKey,
      dnaOrAa,
      selectedGene,
      selectedProtein,
    });
  };
  const onGroupKeyChange = (groupKey) =>
    changeGrouping(groupKey, pending.dnaOrAa);
  const onDnaOrAaChange = (dnaOrAa) =>
    changeGrouping(pending.groupKey, dnaOrAa);

  const getDefaultGeneResidueCoordinates = (selectedGene) => {
    let residueCoordinates = pending.residueCoordinates;
    if (selectedGene.name === 'All Genes') {
      residueCoordinates = [];
    } else {
      residueCoordinates = [[1, selectedGene.len_aa]];
    }
    return residueCoordinates;
  };

  const getDefaultProteinResidueCoordinates = (selectedProtein) => {
    let residueCoordinates = pending.residueCoordinates;
    if (selectedProtein.name === 'All Proteins') {
      residueCoordinates = [];
    } else {
      residueCoordinates = [[1, selectedProtein.len_aa]];
    }
    return residueCoordinates;
  };

  const getCoordinateMode = (coordinateMode) => {
    let { dnaOrAa, residueCoordinates } = pending;
    // If we switched to a coordinate mode that doesn't support AA mutations,
    // then switch off of it now
    if (
      dnaOrAa === DNA_OR_AA.AA &&
      coordinateMode !== COORDINATE_MODES.COORD_GENE &&
      coordinateMode !== COORDINATE_MODES.COORD_PROTEIN
    ) {
      dnaOrAa = DNA_OR_AA.DNA;
    }

    if (coordinateMode === COORDINATE_MODES.COORD_GENE) {
      residueCoordinates = getDefaultGeneResidueCoordinates(
        pending.selectedGene
      );
    } else if (coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
      residueCoordinates = getDefaultProteinResidueCoordinates(
        pending.selectedProtein
      );
    }
    return { dnaOrAa, coordinateMode, residueCoordinates };
  };
  const updateCoordinateMode = (_coordinateMode) => {
    const { dnaOrAa, coordinateMode, residueCoordinates } =
      getCoordinateMode(_coordinateMode);
    setPending({
      ...pending,
      dnaOrAa,
      coordinateMode,
      residueCoordinates,
      validResidueCoordinates: true,
    });
  };

  const updateSelectedGene = (selectedGene) => {
    let { dnaOrAa, coordinateMode, residueCoordinates } = pending;
    if (pending.coordinateMode !== COORDINATE_MODES.COORD_GENE) {
      ({ dnaOrAa, coordinateMode, residueCoordinates } = getCoordinateMode(
        COORDINATE_MODES.COORD_GENE
      ));
    }
    selectedGene = getGene(selectedGene);
    // If we selected a new gene, then update the residue coordinates
    if (selectedGene.name !== pending.selectedGene.name) {
      residueCoordinates = getDefaultGeneResidueCoordinates(selectedGene);
    }
    setPending({
      ...pending,
      dnaOrAa,
      coordinateMode,
      selectedGene,
      residueCoordinates,
      validResidueCoordinates: true,
    });
  };

  const updateSelectedProtein = (selectedProtein) => {
    let { dnaOrAa, coordinateMode, residueCoordinates } = pending;
    if (pending.coordinateMode !== COORDINATE_MODES.COORD_PROTEIN) {
      ({ dnaOrAa, coordinateMode, residueCoordinates } = getCoordinateMode(
        COORDINATE_MODES.COORD_PROTEIN
      ));
    }
    selectedProtein = getProtein(selectedProtein);
    // If we selected a new protein, then update the residue coordinates
    if (selectedProtein.name !== pending.selectedProtein.name) {
      residueCoordinates = getDefaultProteinResidueCoordinates(selectedProtein);
    }
    setPending({
      ...pending,
      dnaOrAa,
      coordinateMode,
      selectedProtein,
      residueCoordinates,
      validResidueCoordinates: true,
    });
  };

  const updateResidueCoordinates = (residueCoordinates) => {
    setPending({
      ...pending,
      residueCoordinates,
      validResidueCoordinates: true,
    });
  };
  const updateValidResidueCoordinates = (valid) => {
    setPending({
      ...pending,
      validResidueCoordinates: valid,
    });
  };

  const updateSelectedPrimers = (selectedPrimers) => {
    let { dnaOrAa, coordinateMode } = pending;
    if (coordinateMode !== COORDINATE_MODES.COORD_PRIMER) {
      ({ dnaOrAa, coordinateMode } = getCoordinateMode(
        COORDINATE_MODES.COORD_PRIMER
      ));
    }
    setPending({
      ...pending,
      dnaOrAa,
      coordinateMode,
      selectedPrimers,
    });
  };

  const updateCustomCoordinates = (customCoordinates) => {
    let { dnaOrAa, coordinateMode } = pending;
    if (coordinateMode !== COORDINATE_MODES.COORD_CUSTOM) {
      ({ dnaOrAa, coordinateMode } = getCoordinateMode(
        COORDINATE_MODES.COORD_CUSTOM
      ));
    }
    setPending({
      ...pending,
      dnaOrAa,
      coordinateMode,
      customCoordinates,
      validCustomCoordinates: true,
    });
  };
  const updateValidCustomCoordinates = (valid) => {
    setPending({
      ...pending,
      validCustomCoordinates: valid,
    });
  };

  const updateCustomSequences = (customSequences) => {
    let { dnaOrAa, coordinateMode } = pending;
    if (coordinateMode !== COORDINATE_MODES.COORD_SEQUENCE) {
      ({ dnaOrAa, coordinateMode } = getCoordinateMode(
        COORDINATE_MODES.COORD_SEQUENCE
      ));
    }
    setPending({
      ...pending,
      dnaOrAa,
      coordinateMode,
      customSequences,
      validCustomSequences: true,
    });
  };
  const updateValidCustomSequences = (valid) => {
    setPending({
      ...pending,
      validCustomSequences: valid,
    });
  };

  const updateSelectedLocationNodes = (selectedLocationNodes) => {
    setPending({
      ...pending,
      selectedLocationNodes,
      // Clear metadata fields
      selectedMetadataFields: {},
    });
  };

  const updateDateRange = (startDate, endDate) => {
    setPending({
      ...pending,
      startDate,
      endDate,
    });
  };

  const updateSubmDateRange = (submStartDate, submEndDate) => {
    setPending({
      ...pending,
      submStartDate,
      submEndDate,
    });
  };

  const updateSelectedMetadataFields = (field, options) => {
    const { selectedMetadataFields } = pending;
    selectedMetadataFields[field] = options;
    setPending({
      ...pending,
      selectedMetadataFields,
    });
  };
  // const updateAgeRange = (ageRange) => {
  //   setPending({
  //     ...pending,
  //     ageRange
  //   });
  // };

  // When the component first mounts (i.e., when the modal is first clicked on)
  // Then:
  //   * Reset the location date tree state to what's currently selected
  //     in the config store
  //   * Fetch metadata fields
  useEffect(() => {
    locationDataStore.setSelectedNodes(configStore.selectedLocationNodes);
    metadataStore.fetchMetadataFields();
  }, []);

  const applyChanges = () => {
    sentRequest.current = true;
    configStore.applyPendingChanges(pending);
  };

  const applyDefault = () => {
    setPending({
      ...pending,
      ...configStore.initialValues,
    });
    // Have to manually trigger this to update the tree
    locationDataStore.setSelectedNodes(
      configStore.initialValues.selectedLocationNodes
    );
  };

  // Make sure everything is in order, before allowing the button to be clicked
  let invalid = false;
  let invalidReason = 'Please fix errors'; // Default message
  if (
    pending.coordinateMode === COORDINATE_MODES.COORD_CUSTOM &&
    !pending.validCustomCoordinates
  ) {
    invalid = true;
    invalidReason = 'Error in custom coordinates';
  } else if (
    pending.coordinateMode === COORDINATE_MODES.COORD_SEQUENCE &&
    !pending.validCustomSequences
  ) {
    invalid = true;
    invalidReason = 'Error in custom sequences';
  } else if (
    (pending.coordinateMode === COORDINATE_MODES.COORD_GENE ||
      pending.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) &&
    !pending.validResidueCoordinates
  ) {
    invalid = true;
    invalidReason = 'Error in residue coordinates';
  } else if (ISOToInt(pending.startDate) > ISOToInt(pending.endDate)) {
    invalid = true;
    invalidReason = 'Start date cannot be before end date';
  } else if (pending.submStartDate !== '' || pending.submStartDate !== '') {
    let submStartDate =
      pending.submStartDate === '' ? MIN_DATE : pending.submStartDate;
    let submEndDate =
      pending.submEndDate === ''
        ? intToISO(new Date().getTime())
        : pending.submEndDate;

    if (ISOToInt(submStartDate) > ISOToInt(submEndDate)) {
      invalid = true;
      invalidReason = 'Start date cannot before end date';
    }
  } else if (pending.selectedLocationNodes.length === 0) {
    invalid = true;
    invalidReason = 'No locations selected';
  } else if (
    pending.coordinateMode === COORDINATE_MODES.COORD_PRIMER &&
    pending.selectedPrimers.length === 0
  ) {
    invalid = true;
    invalidReason = 'No primers selected';
  }
  // When our request goes through, close the modal
  useEffect(() => {
    if (
      sentRequest.current &&
      (UIStore.caseDataState === ASYNC_STATES.SUCCEEDED ||
        UIStore.caseDataState === ASYNC_STATES.FAILED)
    ) {
      // This ref should be unset as when the modal closes it gets wiped from
      // the DOM along with the ref. But just do this in case...
      sentRequest.current = false;
      onRequestClose();
    }
  }, [UIStore.caseDataState]);

  return (
    <Wrapper>
      <Overlay visible={UIStore.caseDataState === ASYNC_STATES.STARTED}>
        <ProgressContainer>
          <LoadingSpinner size={'3rem'} color={'#026cb6'} />
          <ProgressText>Fetching data...</ProgressText>
        </ProgressContainer>
      </Overlay>
      <HeaderContainer>
        <HeaderRow>
          <TitleContainer>
            <div className="title">
              <h2>Select Sequences</h2>
            </div>
          </TitleContainer>
          <div style={{ flexGrow: 1 }} />
          <HeaderButtons>
            {invalid && <InvalidText>Error: {invalidReason}</InvalidText>}
            <CancelButton onClick={onRequestClose}>Cancel</CancelButton>
            <CancelButton onClick={applyDefault}>Reset to Default</CancelButton>
            <ApplyButton
              disabled={invalid}
              invalid={invalid}
              onClick={applyChanges}
            >
              Apply Changes
            </ApplyButton>
          </HeaderButtons>
        </HeaderRow>
      </HeaderContainer>
      <Content>
        <Column minWidth={300} collapseRow={'1'} collapseCol={'1'}>
          <GroupBySelect
            {...pending}
            onGroupKeyChange={onGroupKeyChange}
            onDnaOrAaChange={onDnaOrAaChange}
          />
          <CoordinateSelect
            {...pending}
            updateCoordinateMode={updateCoordinateMode}
            updateSelectedGene={updateSelectedGene}
            updateSelectedProtein={updateSelectedProtein}
            updateResidueCoordinates={updateResidueCoordinates}
            updateValidResidueCoordinates={updateValidResidueCoordinates}
            updateSelectedPrimers={updateSelectedPrimers}
            updateCustomCoordinates={updateCustomCoordinates}
            updateValidCustomCoordinates={updateValidCustomCoordinates}
            updateCustomSequences={updateCustomSequences}
            updateValidCustomSequences={updateValidCustomSequences}
          />
        </Column>
        <Column minWidth={300} collapseRow={'1/3'} collapseCol={'2'}>
          <LocationSelect
            {...pending}
            updateSelectedLocationNodes={updateSelectedLocationNodes}
          />
        </Column>
        <Column minWidth={300} collapseRow={'2'} collapseCol={'1'}>
          <DateSelect
            title={'Date Range'}
            hintText={`<p>Filter for sequences that were collected between the two dates.<br/>Dates are start- and end-inclusive, i.e., [start, end]</p>`}
            startDate={pending.startDate}
            endDate={pending.endDate}
            updateDateRange={updateDateRange}
          />
          <DateSelect
            title={'Submission Date Range'}
            hintText={`<p>Filter for sequences that were <i>submitted</i> to the data repository between the two dates.<br/>Note that submission date does not necessarily equal the date the sequences were published onto the repository (there may be a lag between submission and approval).<br/> Dates are start- and end-inclusive, i.e., [start, end].</p>`}
            startDate={pending.submStartDate}
            endDate={pending.submEndDate}
            updateDateRange={updateSubmDateRange}
          />
          <MetaFieldSelect
            {...pending}
            updateSelectedMetadataFields={updateSelectedMetadataFields}
          />
        </Column>
      </Content>
    </Wrapper>
  );
});

SelectSequencesContent.propTypes = {
  onRequestClose: PropTypes.func,
};
SelectSequencesContent.defaultProps = {
  onRequestClose: NOOP,
};

const SelectSequencesModal = ({ isOpen, onAfterOpen, onRequestClose }) => {
  return (
    <Modal
      isOpen={isOpen}
      onAfterOpen={onAfterOpen}
      onRequestClose={onRequestClose}
      style={{
        overlay: {
          zIndex: 2,
        },
        content: {
          top: '50%',
          left: '50%',
          right: 'auto',
          bottom: 'auto',
          marginRight: '-50%',
          transform: 'translate(-50%, -50%)',
          maxWidth: '100vw',
          zIndex: 3,
          padding: '20px',
          paddingTop: '50px',
          paddingBottom: '0px',
        },
      }}
      contentLabel="Select Sequences"
    >
      <SelectSequencesContent onRequestClose={onRequestClose} />
    </Modal>
  );
};

SelectSequencesModal.propTypes = {
  isOpen: PropTypes.bool.isRequired,
  onAfterOpen: PropTypes.func,
  onRequestClose: PropTypes.func.isRequired,
};
SelectSequencesModal.defaultProps = {
  onAfterOpen: NOOP,
};

export default SelectSequencesModal;
