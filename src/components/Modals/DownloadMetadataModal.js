import React, { useEffect, useRef, useState } from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import PropTypes from 'prop-types';

import {
  ASYNC_STATES,
  SNV_FORMAT,
  GEO_LEVELS,
} from '../../constants/defs.json';

import {
  metadataFields,
  metadataFieldNiceNameMap,
} from '../../constants/metadata';
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
} from './Modal.styles';
import {
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
} from './DownloadMetadataModal.styles';

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

const DownloadMetadataContent = observer(({ onRequestClose }) => {
  const { UIStore, dataStore } = useStores();
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
    snvFormat: SNV_FORMAT.POS_REF_ALT,
  });

  const confirmDownload = () => {
    sentRequest.current = true;
    dataStore.downloadSelectedSequenceMetadata(state);
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

  const handleSnvFormatChange = (event) => {
    setState({
      ...state,
      snvFormat: event.target.value,
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

  let invalid = false;

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

  return (
    <Wrapper>
      <Overlay visible={UIStore.downloadState !== ASYNC_STATES.SUCCEEDED}>
        <ProgressContainer>
          <LoadingSpinner size={'3rem'} color={'#026cb6'} />
          <ProgressText>Fetching data...</ProgressText>
        </ProgressContainer>
      </Overlay>
      <HeaderContainer>
        <HeaderRow>
          <TitleContainer>
            <div className="title">
              <h2>Download Sequence Metadata</h2>
            </div>
          </TitleContainer>
          <div style={{ flexGrow: 1 }} />
          <HeaderButtons>
            <CancelButton onClick={onRequestClose}>Cancel</CancelButton>
            <ApplyButton
              disabled={invalid}
              invalid={invalid}
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
            Downloading sequence metadata for{' '}
            <b>{dataStore.numSequencesAfterAllFiltering}</b> selected sequences
          </Info>
        </Row>
        <Row>
          <RadioForm>
            <FormTitle>SNV Format:</FormTitle>
            <Radio>
              <input
                type="radio"
                value={SNV_FORMAT.POS_REF_ALT}
                checked={state.snvFormat === SNV_FORMAT.POS_REF_ALT}
                onChange={handleSnvFormatChange}
              />
              &lt;Position&gt;|&lt;Reference&gt;|&lt;Alternate&gt; (i.e.,
              &quot;S|614|D|G&quot;)
            </Radio>
            <Radio>
              <input
                type="radio"
                value={SNV_FORMAT.REF_POS_ALT}
                checked={state.snvFormat === SNV_FORMAT.REF_POS_ALT}
                onChange={handleSnvFormatChange}
              />
              &lt;Reference&gt;&lt;Position&gt;&lt;Alternate&gt; (i.e.,
              &quot;S:D614G&quot;)
            </Radio>
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
          <CheckboxForm>
            <FormTitle>SNVs</FormTitle>
            <Checkbox>
              <input
                name="dna"
                type="checkbox"
                checked={state.selectedFields['dna']}
                onChange={handleFieldSelect}
              />
              NT
            </Checkbox>
            <Checkbox>
              <input
                name="gene_aa"
                type="checkbox"
                checked={state.selectedFields['gene_aa']}
                onChange={handleFieldSelect}
              />
              AA (Gene)
            </Checkbox>
            <Checkbox>
              <input
                name="protein_aa"
                type="checkbox"
                checked={state.selectedFields['protein_aa']}
                onChange={handleFieldSelect}
              />
              AA (Protein)
            </Checkbox>
          </CheckboxForm>
        </Row>
      </Content>
    </Wrapper>
  );
});

const DownloadMetadataModal = ({ isOpen, onAfterOpen, onRequestClose }) => {
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
      contentLabel="Download Sequence Metadata"
    >
      <DownloadMetadataContent onRequestClose={onRequestClose} />
    </Modal>
  );
};

DownloadMetadataModal.propTypes = {
  isOpen: PropTypes.bool.isRequired,
  onAfterOpen: PropTypes.func,
  onRequestClose: PropTypes.func.isRequired,
};
DownloadMetadataModal.defaultProps = {
  onAfterOpen: NOOP,
};

export default DownloadMetadataModal;
