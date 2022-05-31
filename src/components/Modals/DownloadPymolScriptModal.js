import React, { useState } from 'react';
import PropTypes from 'prop-types';

import { PYMOL_SCRIPT_TYPES } from '../../constants/defs';

import Modal from 'react-modal';
// import ReactTooltip from 'react-tooltip';
// import QuestionButton from '../Buttons/QuestionButton';

import {
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
  // SelectInput,
  TextInput,
  CheckboxInput,
  ApplyButton,
} from './Modal.styles';
import {} from './DownloadPymolScriptModal.styles';

Modal.setAppElement('#app');
const NOOP = () => {};

const DownloadPymolScriptContent = ({ onRequestClose, onConfirm }) => {
  const [state, setState] = useState({
    scriptType: PYMOL_SCRIPT_TYPES.SCRIPT,
    selectIndividualMutations: true,
    selectAllMutations: true,
    includeDomains: true,
    baseColor: '#FFFFFF',
    ignoreColor: '#AAAAAA',
  });

  // useEffect(() => {
  //   ReactTooltip.rebuild();
  // }, []);

  const confirmDownload = () => {
    onConfirm(state);
  };

  const onChangeScriptType = (event) => {
    setState({
      ...state,
      scriptType: event.target.value,
    });
  };

  const toggleSelectIndividualMutations = (event) => {
    setState({
      ...state,
      selectIndividualMutations: event.target.checked,
    });
  };

  const toggleSelectAllMutations = (event) => {
    setState({
      ...state,
      selectAllMutations: event.target.checked,
    });
  };

  const toggleIncludeDomains = (event) => {
    setState({
      ...state,
      includeDomains: event.target.checked,
    });
  };

  const onChangeBaseColor = (event) => {
    setState({
      ...state,
      baseColor: event.target.value,
    });
  };
  const onChangeIgnoreColor = (event) => {
    setState({
      ...state,
      ignoreColor: event.target.value,
    });
  };

  let invalid = false;
  let invalidReason = '';

  return (
    <Wrapper width={600} height={500}>
      <HeaderContainer>
        <HeaderRow>
          <TitleContainer>
            <div className="title">
              <h2>Download PyMOL Script</h2>
            </div>
          </TitleContainer>
          <div style={{ flexGrow: 1 }} />
          <HeaderButtons>
            {invalid && <InvalidText>Error: {invalidReason}</InvalidText>}
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
            Select &quot;Script&quot; to download a complete Python script,
            which will 1) load the model, 2) select mutations, and 3) color
            mutations by frequency
          </Info>
          <Info>
            Select &quot;Commands&quot; to download a list of PyMOL commands
            that executes the same functions as the script. Select this option
            if you wish to modify certain operations of the script, i.e., adapt
            Spike mutations onto another model of Spike bound to an antibody
          </Info>
        </Row>
        <Row>
          <CheckboxInput>
            <input
              type="radio"
              value={PYMOL_SCRIPT_TYPES.SCRIPT}
              checked={state.scriptType === PYMOL_SCRIPT_TYPES.SCRIPT}
              onChange={onChangeScriptType}
            />
            Script
          </CheckboxInput>
          <CheckboxInput>
            <input
              type="radio"
              value={PYMOL_SCRIPT_TYPES.COMMANDS}
              checked={state.scriptType === PYMOL_SCRIPT_TYPES.COMMANDS}
              onChange={onChangeScriptType}
            />
            Commands
          </CheckboxInput>
        </Row>
        {state.scriptType === PYMOL_SCRIPT_TYPES.SCRIPT && (
          <Row>
            <Info>
              To load this script into PyMOL, open the downloaded .py file with
              PyMOL. (Windows: Open With → Choose another app → Select PyMOL,
              Mac OSX: Open With → Other... → Select PyMOL.app)
            </Info>
            <Info>
              Or, you can first open the PyMOL application and then load this
              script with <code>load path/to/script.py</code>
            </Info>
          </Row>
        )}
        {state.scriptType === PYMOL_SCRIPT_TYPES.COMMANDS && (
          <Row>
            <Info>
              To run these commands, paste the commands from the downloaded file
              into the PyMOL console
            </Info>
          </Row>
        )}
        <Row>
          <CheckboxInput>
            <input
              type="checkbox"
              name="download-pymol-script-select-individual-mutations"
              checked={state.selectIndividualMutations}
              onChange={toggleSelectIndividualMutations}
            />
            Create a selection for each individual mutation
          </CheckboxInput>
        </Row>
        <Row>
          <CheckboxInput>
            <input
              type="checkbox"
              name="download-pymol-script-select-all-mutations"
              checked={state.selectAllMutations}
              onChange={toggleSelectAllMutations}
            />
            Create a selection for all mutations
          </CheckboxInput>
        </Row>
        <Row>
          <CheckboxInput>
            <input
              type="checkbox"
              name="download-pymol-script-include-domains"
              checked={state.includeDomains}
              onChange={toggleIncludeDomains}
            />
            Create a selection for protein domains
          </CheckboxInput>
        </Row>
        <Row>
          <TextInput>
            Heatmap Base Color
            <input
              type="color"
              placeholder="#RRGGBB"
              value={state.baseColor}
              onChange={onChangeBaseColor}
            />
            <span style={{ marginLeft: 5 }}>{state.baseColor}</span>
          </TextInput>
        </Row>
        <Row>
          <TextInput>
            Non-heatmap Base Color
            <br />
            (for chains where heatmap is not applied)
            <input
              type="color"
              placeholder="#RRGGBB"
              value={state.ignoreColor}
              onChange={onChangeIgnoreColor}
            />
            <span style={{ marginLeft: 5 }}>{state.ignoreColor}</span>
          </TextInput>
        </Row>
      </Content>
    </Wrapper>
  );
};
DownloadPymolScriptContent.propTypes = {
  onRequestClose: PropTypes.func.isRequired,
  onConfirm: PropTypes.func.isRequired,
};

const DownloadPymolScriptModal = ({
  isOpen,
  onAfterOpen,
  onRequestClose,
  onConfirm,
}) => {
  const closeDownloadModal = () => {
    onRequestClose();
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
      contentLabel="Download PyMOL Script"
    >
      <DownloadPymolScriptContent
        onRequestClose={closeDownloadModal}
        onConfirm={onConfirm}
      />
    </Modal>
  );
};

DownloadPymolScriptModal.propTypes = {
  isOpen: PropTypes.bool.isRequired,
  onAfterOpen: PropTypes.func,
  onRequestClose: PropTypes.func.isRequired,
  onConfirm: PropTypes.func.isRequired,
};
DownloadPymolScriptModal.defaultProps = {
  onAfterOpen: NOOP,
};

export default DownloadPymolScriptModal;
