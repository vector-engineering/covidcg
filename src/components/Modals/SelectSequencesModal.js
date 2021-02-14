import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';

import Button from '../Buttons/Button';
import Modal from 'react-modal';

import LocationSelect from '../Selection/LocationSelect';
import GroupBySelect from '../Selection/GroupBySelect';
import CoordinateSelect from '../Selection/CoordinateSelect';
import DateSelect from '../Selection/DateSelect';
import MetaFieldSelect from '../Selection/MetaFieldSelect';
import FilterDataIntoOther from '../Selection/FilterDataIntoOther';

const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;
  justify-content: flex-start;

  width: calc(100vw - 100px);
  height: calc(100vh - 100px);
`;

const Footer = styled.div`
  height: 40px;

  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-end;

  border-top: 1px solid #ccc;

  padding: 5px 10px;
  margin-left: -20px;
  margin-right: -20px;
`;

const Header = styled.div`
  display: flex;
  flex-direction: row;
  align-items: flex-start;

  padding: 10px;
  margin-right: 15px;
  margin-bottom: 10px;
  border: 1px solid #ccc;

  .title {
    h2 {
      margin-bottom: 0px;
      margin-top: 0px;
    }
  }
  .spacer {
    flex-grow: 1;
  }
  .close-button {
  }
`;
const Content = styled.div`
  height: 100%;

  font-size: 1em;
  font-weight: normal;

  display: flex;
  flex-direction: row;
  align-items: flex-start;
  flex-wrap: wrap;
`;

const Column = styled.div`
  width: ${({ width }) => width}px;
  height: 100%;

  display: flex;
  flex-direction: column;
  align-items: stretch;
  overflow-y: scroll;
  ${({ grow }) => (grow ? 'flex-grow: 1;' : '')}
`;
Column.defaultProps = {
  width: 300,
  grow: false,
};

const CancelButton = styled(Button)`
  background-color: #ddd;
  color: #000;
  background-image: none;
  font-size: 1rem;
  font-weight: normal;
  margin-right: 10px;

  &:hover,
  &:active {
    background-color: #eee;
  }
`;

const ApplyButton = styled(Button)`
  font-size: 1rem;
  font-weight: normal;
`;

Modal.setAppElement('#app');
const NOOP = () => {};

const SelectSequencesContent = ({ onRequestClose }) => {
  const applyChanges = () => {};

  return (
    <Wrapper>
      <Content>
        <Column width={300}>
          <Header>
            <div className="title">
              <h2>Select Sequences</h2>
            </div>
          </Header>
          <GroupBySelect />
          <CoordinateSelect />
        </Column>
        <Column width={300}>
          <LocationSelect />
        </Column>
        <Column width={300} grow={true}>
          <DateSelect />
          <MetaFieldSelect />
          <FilterDataIntoOther />
        </Column>
      </Content>
      <Footer>
        <CancelButton onClick={onRequestClose}>Cancel</CancelButton>
        <ApplyButton onClick={applyChanges}>Apply Changes</ApplyButton>
      </Footer>
    </Wrapper>
  );
};

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
          padding: '20px 20px 0px 20px',
          zIndex: 3,
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
