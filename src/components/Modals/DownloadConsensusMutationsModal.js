import React, { useEffect, useState, useRef } from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { config } from '../../config';

import { ASYNC_STATES } from '../../constants/defs.json';
import { getReferenceNames, getReferences } from '../../utils/reference';

import Modal from 'react-modal';
import ReactTooltip from 'react-tooltip';
import QuestionButton from '../Buttons/QuestionButton';

import LoadingSpinner from '../Common/LoadingSpinner';

import {
  Wrapper,
  Overlay,
  Content,
  Row,
  ProgressContainer,
  ProgressText,
  TitleContainer,
  HeaderContainer,
  HeaderRow,
  HeaderButtons,
  CancelButton,
  InvalidText,
  ApplyButton,
  SelectInput,
  TextInput,
} from './Modal.styles';

Modal.setAppElement('#app');
const NOOP = () => {};

const DownloadConsensusMutationsContent = observer(({ onRequestClose }) => {
  const { UIStore, groupDataStore } = useStores();
  const sentRequest = useRef();
  const [state, setState] = useState({
    group: Object.keys(config.group_cols)[0],
    mutationType: 'gene_aa',
    consensusThreshold: 0.9,
    selectedReference: groupDataStore.reportActiveReference,
  });

  useEffect(() => {
    ReactTooltip.rebuild();
  }, []);

  const confirmDownload = () => {
    sentRequest.current = true;
    groupDataStore.downloadGroupMutationFrequencyData({
      group: state.group,
      mutationType: state.mutationType,
      consensusThreshold: parseFloat(state.consensusThreshold),
      selectedReference: state.selectedReference,
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
  let invalidReason = '';

  if (UIStore.downloadState === ASYNC_STATES.FAILED) {
    invalid = true;
    invalidReason = 'Download Failed';
  }

  const onChangeGroup = (e) => {
    setState({
      ...state,
      group: e.target.value,
    });
  };

  const onChangeReference = (e) => {
    setState({
      ...state,
      selectedReference: e.target.value,
    });
  };

  const onChangeMutationType = (e) => {
    setState({
      ...state,
      mutationType: e.target.value,
    });
  };

  const onChangeConsensusThreshold = (e) => {
    setState({
      ...state,
      consensusThreshold: e.target.value,
    });
  };

  const groupOptions = [];
  Object.keys(config.group_cols).forEach((group) => {
    groupOptions.push(
      <option key={`group-${group}`} value={group}>
        {config.group_cols[group]['title']}
      </option>
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

  return (
    <Wrapper>
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
              <h2>Download Consensus Mutations</h2>
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
          <SelectInput>
            Phylogeny Definition{' '}
            <select onChange={onChangeGroup} value={state.group}>
              {groupOptions}
            </select>
          </SelectInput>
        </Row>
        <Row>
          <SelectInput>
            Reference{' '}
            <select
              onChange={onChangeReference}
              value={state.selectedReference}
              style={{ width: '100%' }}
            >
              {referenceOptionItems}
            </select>
          </SelectInput>
        </Row>
        <Row>
          <SelectInput>
            Mutation Type{' '}
            <select onChange={onChangeMutationType} value={state.mutationType}>
              <option value="dna">NT</option>
              <option value="gene_aa">AA (Gene)</option>
              <option value="protein_aa">AA (Protein)</option>
            </select>
          </SelectInput>
        </Row>
        <Row>
          <TextInput>
            Consensus Threshold
            <QuestionButton
              data-tip={`<p>Mutations with frequencies below this percentage are excluded from the download. For example, if a mutation is only present in 80% of sequences in a given lineage, and the Consensus Threshold is set to 0.9, then the mutation will be excluded.</p>`}
              data-html="true"
              data-place="right"
              data-for="main-tooltip"
            />
            <input
              type="number"
              min={0}
              max={1}
              step={0.01}
              value={state.consensusThreshold}
              onChange={onChangeConsensusThreshold}
            />
          </TextInput>
        </Row>
      </Content>
    </Wrapper>
  );
});

const DownloadConsensusMutationsModal = ({
  isOpen,
  onAfterOpen,
  onRequestClose,
}) => {
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
      contentLabel="Download Sequence Metadata"
    >
      <DownloadConsensusMutationsContent onRequestClose={closeDownloadModal} />
    </Modal>
  );
};

DownloadConsensusMutationsModal.propTypes = {
  isOpen: PropTypes.bool.isRequired,
  onAfterOpen: PropTypes.func,
  onRequestClose: PropTypes.func.isRequired,
};
DownloadConsensusMutationsModal.defaultProps = {
  onAfterOpen: NOOP,
};

export default DownloadConsensusMutationsModal;
