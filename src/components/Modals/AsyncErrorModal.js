import React from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import PropTypes from 'prop-types';

import Modal from 'react-modal';
import ExternalLink from '../Common/ExternalLink';

import {
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
  RefreshButton,
} from './AsyncErrorModal.styles';

Modal.setAppElement('#app');
const NOOP = () => {};

const AsyncErrorContent = observer(({ onRequestClose }) => {
  // const { UIStore, dataStore } = useStores();
  let invalid = false;

  const refreshPage = () => {
    window.location.reload();
  };

  return (
    <Wrapper>
      <HeaderContainer>
        <HeaderRow>
          <TitleContainer>
            <div className="title">
              <h2>Error Fetching Data</h2>
            </div>
          </TitleContainer>
          <div style={{ flexGrow: 1 }} />
          <HeaderButtons>
            <CancelButton onClick={onRequestClose}>Close</CancelButton>
            <RefreshButton
              disabled={invalid}
              invalid={invalid}
              onClick={refreshPage}
            >
              Refresh Page
            </RefreshButton>
          </HeaderButtons>
        </HeaderRow>
      </HeaderContainer>
      <Content>
        <Row>
          <Info>
            There was an error getting data from the server. To try again, first
            close this box by pressing &quot;Cancel&quot; at the top, or by
            clicking anywhere outside the box
          </Info>
        </Row>
        <Row>
          <Info>
            If this error persists, please contact us at{' '}
            <ExternalLink href="mailto:covidcg@broadinstitute.org">
              covidcg@broadinstitute.org
            </ExternalLink>
          </Info>
        </Row>
      </Content>
    </Wrapper>
  );
});

const AsyncErrorModal = ({ isOpen, onAfterOpen, onRequestClose }) => {
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
      contentLabel="Error fetching data"
    >
      <AsyncErrorContent onRequestClose={onRequestClose} />
    </Modal>
  );
};

AsyncErrorModal.propTypes = {
  isOpen: PropTypes.bool.isRequired,
  onAfterOpen: PropTypes.func,
  onRequestClose: PropTypes.func.isRequired,
};
AsyncErrorModal.defaultProps = {
  onAfterOpen: NOOP,
};

export default AsyncErrorModal;
