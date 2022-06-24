import React, { useEffect, useRef, useState } from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import PropTypes from 'prop-types';

import {
  ASYNC_STATES,
  MUTATION_FORMAT,
  GEO_LEVELS,
} from '../../constants/defs.json';

import {
  metadataFields,
  metadataFieldNiceNameMap,
} from '../../constants/metadata';
import { getReferenceNames, getReferences } from '../../utils/reference';
import { config } from '../../config';

import Modal from 'react-modal';

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
  Wrapper,
  Content,
  Row,
  Info,
  RadioForm,
  Radio,
  CheckboxForm,
  FormTitle,
  Checkbox,
  ApplyButton,
  SelectInput,
} from './Modal.styles';

Modal.setAppElement('#app');
const NOOP = () => {};

const initialSelectedMetadataFields = {};
metadataFields.forEach((field) => {
  initialSelectedMetadataFields[field] = true;
});

const initialSelectedGroupings = {};
Object.keys(config.group_cols).forEach((group) => {
  initialSelectedGroupings[group] = true;
});

const initialSelectedLocationFields = {};
Object.values(GEO_LEVELS).forEach((field) => {
  initialSelectedLocationFields[field] = true;
});

const DownloadVariantTableContent = observer(({ onRequestClose }) => {
  const { UIStore, dataStore, configStore } = useStores();
  const sentRequest = useRef();

  const [state, setState] = useState({
    selectedFields: Object.assign(
      {
        dna: true,
        gene_aa: true,
        protein_aa: true,
      },
      initialSelectedMetadataFields,
      initialSelectedGroupings,
      initialSelectedLocationFields
    ),
    selectedReference: configStore.selectedReference,
    mutationFormat: MUTATION_FORMAT.REF_POS_ALT,
  });

  const confirmDownload = () => {
    sentRequest.current = true;
    dataStore.downloadVariantTable(state);
  };

  const handleFieldSelect = (event) => {
    const target = event.target;

    const selectedFields = Object.assign({}, state.selectedFields);
    selectedFields[target.name] = target.checked;

    setState({
      ...state,
      selectedFields,
    });
  };

  const handleMutationFormatChange = (event) => {
    setState({
      ...state,
      mutationFormat: event.target.value,
    });
  };

  const handleReferenceChange = (e) => {
    setState({
      ...state,
      selectedReference: e.target.value,
    });
  };

  // When our request goes through, close the modal
  useEffect(() => {
    if (
      sentRequest.current &&
      UIStore.downloadState === ASYNC_STATES.SUCCEEDED
    ) {
      // This ref should be unset as when the modal closes it gets wiped from
      // the DOM along with the ref. But just do this in case...
      sentRequest.current = false;
      onRequestClose();
    }
  }, [UIStore.downloadState]);

  // Build all of the checkboxes
  const fieldCheckboxes = [];
  metadataFields.forEach((field) => {
    const fieldNiceName = metadataFieldNiceNameMap[field];
    fieldCheckboxes.push(
      <Checkbox key={`metadata-download-metafield-${field}`}>
        <input
          name={field}
          type="checkbox"
          checked={state.selectedFields[field]}
          onChange={handleFieldSelect}
        />
        {fieldNiceName}
      </Checkbox>
    );
  });

  const groupingCheckboxes = [];
  Object.keys(config.group_cols).forEach((group) => {
    groupingCheckboxes.push(
      <Checkbox key={`metadata-download-grouping-${group}`}>
        <input
          name={group}
          type="checkbox"
          checked={state.selectedFields[group]}
          onChange={handleFieldSelect}
        />
        {config.group_cols[group].title}
      </Checkbox>
    );
  });

  // Location fields
  const locationCheckboxes = [];
  Object.values(GEO_LEVELS).forEach((field) => {
    locationCheckboxes.push(
      <Checkbox key={`metadata-download-location-${field}`}>
        <input
          name={field}
          type="checkbox"
          checked={state.selectedFields[field]}
          onChange={handleFieldSelect}
        />
        {/* Capitalize */}
        {field.charAt(0).toUpperCase() + field.slice(1)}
      </Checkbox>
    );
  });

  const referenceOptionItems = [];
  getReferenceNames().forEach((referenceName) => {
    let name =
      referenceName + ' ' + getReferences()[referenceName]['description'];

    referenceOptionItems.push(
      <option key={`ref-option-${referenceName}`} value={referenceName}>
        {name}
      </option>
    );
  });

  let invalid = false;
  let invalidReason = '';

  if (UIStore.downloadState === ASYNC_STATES.FAILED) {
    invalid = true;
    invalidReason = 'Download Failed';
  }

  return (
    <Wrapper width={600} height={400}>
      <Overlay visible={UIStore.downloadState === ASYNC_STATES.STARTED}>
        <ProgressContainer>
          <LoadingSpinner size={'3rem'} color={'#026cb6'} />
          <ProgressText>Fetching data...</ProgressText>
        </ProgressContainer>
      </Overlay>
      <HeaderContainer>
        <HeaderRow>
          <TitleContainer>
            <div className="title">
              <h2>Download Variant Table</h2>
            </div>
          </TitleContainer>
          <div style={{ flexGrow: 1 }} />
          <HeaderButtons>
            {invalid && <InvalidText>Error: {invalidReason}</InvalidText>}
            <CancelButton onClick={onRequestClose}>Cancel</CancelButton>
            <ApplyButton
              disabled={UIStore.downloadState === ASYNC_STATES.STARTED}
              invalid={UIStore.downloadState === ASYNC_STATES.STARTED}
              onClick={confirmDownload}
            >
              Download
            </ApplyButton>
          </HeaderButtons>
        </HeaderRow>
      </HeaderContainer>
      <Content>
        <Row>
          <Info>
            Downloading variant table for{' '}
            <b>{dataStore.numSequencesAfterAllFiltering}</b> selected sequences
          </Info>
        </Row>
        <Row>
          <RadioForm>
            <FormTitle>Mutation Format:</FormTitle>
            <Radio>
              <input
                type="radio"
                value={MUTATION_FORMAT.POS_REF_ALT}
                checked={state.mutationFormat === MUTATION_FORMAT.POS_REF_ALT}
                onChange={handleMutationFormatChange}
              />
              &lt;Position&gt;|&lt;Reference&gt;|&lt;Alternate&gt; (i.e.,
              &quot;S|614|D|G&quot;)
            </Radio>
            <Radio>
              <input
                type="radio"
                value={MUTATION_FORMAT.REF_POS_ALT}
                checked={state.mutationFormat === MUTATION_FORMAT.REF_POS_ALT}
                onChange={handleMutationFormatChange}
              />
              &lt;Reference&gt;&lt;Position&gt;&lt;Alternate&gt; (i.e.,
              &quot;S:D614G&quot;)
            </Radio>
          </RadioForm>
        </Row>
        <Row>
          <RadioForm>
            <FormTitle>Reference:</FormTitle>
            <SelectInput>
              <select
                onChange={handleReferenceChange}
                value={state.selectedReference}
                style={{ width: '100%' }}
              >
                {referenceOptionItems}
              </select>
            </SelectInput>
          </RadioForm>
        </Row>
        <Row>
          <CheckboxForm>
            <FormTitle>Metadata Fields</FormTitle>
            {fieldCheckboxes}
          </CheckboxForm>
          <CheckboxForm>
            <FormTitle>Phylogeny</FormTitle>
            {groupingCheckboxes}
          </CheckboxForm>
          <CheckboxForm>
            <FormTitle>Geography</FormTitle>
            {locationCheckboxes}
          </CheckboxForm>
        </Row>
      </Content>
    </Wrapper>
  );
});

const DownloadVariantTableModal = ({ isOpen, onAfterOpen, onRequestClose }) => {
  const { UIStore } = useStores();

  const closeDownloadModal = () => {
    onRequestClose();
    // Before we close the modal, clear the download state
    if (UIStore.downloadState !== ASYNC_STATES.STARTED) {
      UIStore.onDownloadFinished();
    }
  };

  return (
    <Modal
      isOpen={isOpen}
      onAfterOpen={onAfterOpen}
      onRequestClose={closeDownloadModal}
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
      contentLabel="Download Variant Table"
    >
      <DownloadVariantTableContent onRequestClose={closeDownloadModal} />
    </Modal>
  );
};

DownloadVariantTableModal.propTypes = {
  isOpen: PropTypes.bool.isRequired,
  onAfterOpen: PropTypes.func,
  onRequestClose: PropTypes.func.isRequired,
};
DownloadVariantTableModal.defaultProps = {
  onAfterOpen: NOOP,
};

export default DownloadVariantTableModal;
