import React, { useEffect, useRef } from 'react';
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

  const confirmDownload = () => {
    sentRequest.current = true;
    dataStore.downloadGenomes();
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
              <h2>Download Genomes</h2>
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
            Downloading genomes for{' '}
            <b>{dataStore.numSequencesAfterAllFiltering}</b> selected sequences
          </Info>
        </Row>
        <Row>
          <Info>
            <b>Note</b>: downloading a large number of genomes (i.e., 10,000+)
            may take a few minutes
          </Info>
        </Row>
      </Content>
    </Wrapper>
  );
});

const DownloadGenomesModal = ({ isOpen, onAfterOpen, onRequestClose }) => {
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
      <DownloadGenomesContent onRequestClose={onRequestClose} />
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
