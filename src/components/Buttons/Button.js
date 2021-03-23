import React from 'react';
import PropTypes from 'prop-types';

import { StyledButton } from './Button.styles';

const Button = ({
  children,
  sticky,
  disabled,
  onClick,
  onFocus,
  onBlur,
  ...rest
}) => {
  const handleClick = (e) => {
    if (disabled) {
      return;
    }
    onClick(e);
  };
  const handleFocus = () => onFocus();
  const handleBlur = (e) => onBlur(e);

  return (
    <StyledButton
      type="button"
      onClick={handleClick}
      onFocus={handleFocus}
      onBlur={handleBlur}
      sticky={sticky}
      disabled={disabled}
      {...rest}
    >
      {children}
    </StyledButton>
  );
};

Button.propTypes = {
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node,
  ]),
  sticky: PropTypes.bool,
  disabled: PropTypes.bool,
  onClick: PropTypes.func,
  onFocus: PropTypes.func,
  onBlur: PropTypes.func,
};

Button.defaultProps = {
  children: [],
  sticky: false,
  disabled: false,
  // no-op
  onClick: () => {},
  onFocus: () => {},
  onBlur: () => {},
};

export default Button;
