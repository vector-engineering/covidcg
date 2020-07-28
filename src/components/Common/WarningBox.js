import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';

const WarningContainer = styled.div`
  display: ${({ show }) => (show ? 'flex' : 'none')};
  flex-direction: column;
  align-items: stretch;
  justify-content: flex-start;

  // colors from Bootstrap
  background-color: #fff3cd;
  border: 1px solid #aaa;
  border-radius: 5px;

  margin: 0px 10px;
  margin-bottom: 15px;
  padding: 10px 20px;

  .warning-header {
    display: flex;
    flex-direction: row;
    align-items: center;
    margin-bottom: 5px;
    .warning-title {
      font-size: 1.25em;
    }
    .spacer {
      flex-grow: 1;
    }
  }

  .warning-text {
    margin: 0px;
    font-weight: normal;
  }
`;
WarningContainer.defaultProps = {
  show: true,
};

const WarningBox = ({ title, text, show, onDismiss }) => {
  const onDismissWarning = (event) => {
    event.preventDefault();
    onDismiss();
  };

  return (
    <WarningContainer show={show}>
      <div className="warning-header">
        <span className="warning-title">{title}</span>
        <div className="spacer" />
        <button className="close-button" onClick={onDismissWarning}>
          Dismiss
        </button>
      </div>
      <p className="warning-text">{text}</p>
    </WarningContainer>
  );
};
WarningBox.propTypes = {
  title: PropTypes.string,
  text: PropTypes.string,
  show: PropTypes.bool.isRequired,
  onDismiss: PropTypes.func.isRequired,
};
WarningBox.defaultProps = {
  title: 'WARNING:',
  text: 'Default warning text',
};

export default WarningBox;
