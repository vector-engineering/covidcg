import React, { useEffect, useRef } from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import PropTypes from 'prop-types';
import { COORDINATE_MODES, ASYNC_STATES } from '../../constants/defs.json';

import Modal from 'react-modal';

import LocationSelect from '../Selection/LocationSelect';
import GroupBySelect from '../Selection/GroupBySelect';
import CoordinateSelect from '../Selection/CoordinateSelect';
import DateSelect from '../Selection/DateSelect';
import MetaFieldSelect from '../Selection/MetaFieldSelect';
import FilterDataIntoOther from '../Selection/FilterDataIntoOther';

import LoadingSpinner from '../Common/LoadingSpinner';

import {
  Wrapper,
  Content,
  Column,
  Header,
  Footer,
  CancelButton,
  ApplyButton,
  InvalidText,
  Overlay,
  ProgressContainer,
  ProgressText,
} from './SelectSequencesModal.styles';

Modal.setAppElement('#app');
const NOOP = () => {};

const SelectSequencesContent = observer(({ onRequestClose }) => {
  const { configStore, dataStore, UIStore } = useStores();
  const sentRequest = useRef(false);

  const applyChanges = () => {
    sentRequest.current = true;
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

  // When our request goes through, close the modal
  useEffect(() => {
    if (
      sentRequest.current &&
      UIStore.caseDataState === ASYNC_STATES.SUCCEEDED
    ) {
      // This ref should be unset as when the modal closes it gets wiped from
      // the DOM along with the ref. But just do this in case...
      sentRequest.current = false;
      onRequestClose();
    }
  }, [UIStore.caseDataState]);

  return (
    <Wrapper>
      <Overlay visible={UIStore.caseDataState !== ASYNC_STATES.SUCCEEDED}>
        <ProgressContainer>
          <LoadingSpinner size={'3rem'} color={'#026cb6'} />
          <ProgressText>Fetching data...</ProgressText>
        </ProgressContainer>
      </Overlay>
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
