import React, { useEffect, useRef, useState } from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import PropTypes from 'prop-types';

import { ASYNC_STATES } from '../../constants/defs.json';

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
  CheckboxInput,
} from './Modal.styles';
import {
  Wrapper,
  Content,
  Row,
  Info,
  ApplyButton,
} from './DownloadGenomesModal.styles';

Modal.setAppElement('#app');
const NOOP = () => {};

const DownloadGenomesContent = observer(({ onRequestClose }) => {
  const { UIStore, dataStore } = useStores();
  const sentRequest = useRef();
  const [state, setState] = useState({
    compress: true,
  });

  const handleCompressSelect = (e) => {
    setState({
      compress: e.target.checked,
    });
  };

  const confirmDownload = () => {
    sentRequest.current = true;
    dataStore.downloadGenomes({ compress: state.compress });
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
              <h2>Download Genomes</h2>
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
            Downloading genomes for{' '}
            <b>{dataStore.numSequencesAfterAllFiltering}</b> selected sequences
          </Info>
        </Row>
        <Row>
          <CheckboxInput>
            <input
              name="compress"
              type="checkbox"
              checked={state.compress}
              onChange={handleCompressSelect}
            ></input>
            Compress Output
          </CheckboxInput>
        </Row>
        <Row>
          <Info>
            <b>Note</b>: downloading a large number of genomes (i.e., 10,000+)
            may take a few minutes
          </Info>
          <Info>
            Selecting &quot;Compress Output&quot; will compress the resulting
            FASTA file on our server. This may take much longer but will reduce
            the filesize and bandwidth required.
          </Info>
        </Row>
      </Content>
    </Wrapper>
  );
});

const DownloadGenomesModal = ({ isOpen, onAfterOpen, onRequestClose }) => {
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
      <DownloadGenomesContent onRequestClose={closeDownloadModal} />
    </Modal>
  );
};

DownloadGenomesModal.propTypes = {
  isOpen: PropTypes.bool.isRequired,
  onAfterOpen: PropTypes.func,
  onRequestClose: PropTypes.func.isRequired,
};
DownloadGenomesModal.defaultProps = {
  onAfterOpen: NOOP,
};

export default DownloadGenomesModal;
