import React from 'react';
import PropTypes from 'prop-types';

import {
  WarningContainer,
  WarningHeader,
  WarningText,
} from './WarningBox.styles';

const WarningBox = ({
  children,
  title,
  show,
  onDismiss,
  showDismissButton,
}) => {
  const onDismissWarning = (event) => {
    event.preventDefault();
    onDismiss();
  };

  return (
    <WarningContainer show={show}>
      <WarningHeader>
        <span className="warning-title">{title}</span>
        <div className="spacer" />
        {showDismissButton && (
          <button className="close-button" onClick={onDismissWarning}>
            Dismiss
          </button>
        )}
      </WarningHeader>
      <WarningText>{children}</WarningText>
    </WarningContainer>
  );
};
WarningBox.propTypes = {
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node,
  ]),
  title: PropTypes.string,
  show: PropTypes.bool.isRequired,
  onDismiss: PropTypes.func,
  showDismissButton: PropTypes.bool,
};
WarningBox.defaultProps = {
  title: 'WARNING:',
  children: 'Default warning text',
  onDismiss: () => {},
  showDismissButton: true,
};

export default WarningBox;
