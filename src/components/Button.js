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
`;

const Button = ({ children, onClick }) => {
  const handleClick = () => {
    onClick();
  };

  return <StyledButton onClick={handleClick}>{children}</StyledButton>;
};

Button.propTypes = {
  children: PropTypes.oneOfType([PropTypes.element, PropTypes.string]),
  onClick: PropTypes.func,
};

Button.defaultProps = {
  children: [],
  // no-op
  onClick: () => {},
};

export default Button;
