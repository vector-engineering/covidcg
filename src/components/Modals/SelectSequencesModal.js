import React from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import PropTypes from 'prop-types';
import { COORDINATE_MODES } from '../../constants/defs.json';

import Modal from 'react-modal';

import LocationSelect from '../Selection/LocationSelect';
import GroupBySelect from '../Selection/GroupBySelect';
import CoordinateSelect from '../Selection/CoordinateSelect';
import DateSelect from '../Selection/DateSelect';
import MetaFieldSelect from '../Selection/MetaFieldSelect';
import FilterDataIntoOther from '../Selection/FilterDataIntoOther';

import {
  Wrapper,
  Content,
  Column,
  Header,
  Footer,
  CancelButton,
  ApplyButton,
  InvalidText,
} from './SelectSequencesModal.styles';

Modal.setAppElement('#app');
const NOOP = () => {};

const SelectSequencesContent = observer(({ onRequestClose }) => {
  const { configStore, dataStore } = useStores();

  const applyChanges = () => {
    dataStore.fetchData();
  };

  // Make sure everything is in order, before allowing the button to be clicked
  let invalid = false;
  if (
    configStore.coordinateMode === COORDINATE_MODES.COORD_CUSTOM &&
    !configStore.validCustomCoordinates
  ) {
    invalid = true;
  } else if (
    configStore.coordinateMode === COORDINATE_MODES.COORD_SEQUENCE &&
    !configStore.validCustomSequences
  ) {
    invalid = true;
  } else if (
    (configStore.coordinateMode === COORDINATE_MODES.COORD_GENE ||
      configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) &&
    !configStore.validResidueCoordinates
  ) {
    invalid = true;
  } else if (!configStore.validDateRange) {
    invalid = true;
  }

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
        {invalid && <InvalidText>Please fix errors</InvalidText>}
        <CancelButton onClick={onRequestClose}>Cancel</CancelButton>
        <ApplyButton
          disabled={invalid}
          invalid={invalid}
          onClick={applyChanges}
        >
          Apply Changes
        </ApplyButton>
      </Footer>
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
