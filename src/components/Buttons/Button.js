import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';

const StyledButton = styled.button`
  padding: 5px 10px;
  background-color: #28a745;
  background-image: linear-gradient(-180deg, #34d058, #28a745 90%);
  color: #fff;
  font-size: 12px;
  border: 1px solid rgba(27, 31, 35, 0.2);
  border-radius: 0.25em;
  outline: none;

  &:active {
    background-color: #279f43;
    background-image: none;
    border-color: rgba(27, 31, 35, 0.5);
  }

  ${({ sticky }) =>
    sticky &&
    `
    &:focus {
      background-color: #279f43;
      background-image: none;
      border-color: rgba(27, 31, 35, 0.5);
    }
  `}
`;

const Button = ({ children, sticky, onClick, onFocus, onBlur, ...rest }) => {
  const handleClick = (e) => onClick(e);
  const handleFocus = () => onFocus();
  const handleBlur = (e) => onBlur(e);

  return (
    <StyledButton
      type="button"
      onClick={handleClick}
      onFocus={handleFocus}
      onBlur={handleBlur}
      sticky={sticky}
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
  onClick: PropTypes.func,
  onFocus: PropTypes.func,
  onBlur: PropTypes.func,
};

Button.defaultProps = {
  children: [],
  sticky: false,
  // no-op
  onClick: () => {},
  onFocus: () => {},
  onBlur: () => {},
};

export default Button;
