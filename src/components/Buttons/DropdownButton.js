import React, { useState } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import _ from 'underscore';

import Button from './Button';

const DropdownContainer = styled.div`
  position: relative;

  button {
    padding-right: 18px;
  }

  .caret {
    position: relative;
    margin-left: 5px;
  }

  .caret:after {
    content: '';
    position: absolute;
    left: 1px;
    top: 5px;
    border-top: 5px solid #eeeeee;
    border-left: 5px solid transparent;
    border-right: 5px solid transparent;
  }
`;

const DropdownMenu = styled.div`
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  justify-content: center;

  position: absolute;
  top: 28px;
  right: 0;
  z-index: 10;

  min-width: 10rem;
  padding: 5px 0px;

  background-color: #fff;
  border: 1px solid #aaa;
  border-radius: 3px;

  overflow: hidden;
`;

// DropdownMenu.defaultProps = {
//   direction: 'left',
// };

const DropdownItem = styled.a`
  display: block;
  width: 100%;
  padding: 0.25rem 1.5rem;
  clear: both;
  font-weight: 400;
  color: #212529;
  text-align: inherit;
  white-space: nowrap;
  background-color: #fff;
  border: 0;
  text-decoration: none;

  &:hover {
    background-color: #f8f9fa;
  }
  &:active {
    color: #ffffff;
    background-color: #34d058;
  }
`;

const DropdownButton = ({
  button,
  text,
  options,
  values,
  onSelect,
  ...other
}) => {
  const [state, setState] = useState({
    expanded: false,
    id: _.uniqueId('dropdown_'),
  });

  const handleToggle = () => {
    setState({ ...state, expanded: !state.expanded });
  };

  // const handleFocus = () => {
  //   console.log('focus');
  // };

  const handleBlur = (e) => {
    if (
      e.relatedTarget == null ||
      !e.relatedTarget.classList.contains('dropdown-item')
    ) {
      setState({ ...state, expanded: false });
    }
  };
  // const handleBlur = () => {};

  const handleOption = (e) => {
    onSelect(e.target.dataset.option);
    setState({ ...state, expanded: false });
  };

  // Build option anchor links
  const optionElements = [];
  if (values === null) {
    values = options;
  }
  for (let i = 0; i < options.length; i++) {
    const option = options[i];
    const value = values[i];
    optionElements.push(
      <DropdownItem
        key={option}
        className="dropdown-item"
        href="#"
        onClick={handleOption}
        data-option={value}
      >
        {option}
      </DropdownItem>
    );
  }

  const ButtonElement = button;

  return (
    <DropdownContainer>
      <ButtonElement
        id={state.id}
        data-toggle="dropdown"
        aria-haspopup="true"
        aria-expanded={state.expanded}
        onClick={handleToggle}
        // onFocus={handleFocus}
        onBlur={handleBlur}
        sticky={true}
        {...other}
      >
        {text}
        <span className="caret"></span>
      </ButtonElement>
      <DropdownMenu
        aria-labelledby={state.id}
        style={{
          display: state.expanded ? 'flex' : 'none',
        }}
      >
        {optionElements}
      </DropdownMenu>
    </DropdownContainer>
  );
};

DropdownButton.propTypes = {
  button: PropTypes.elementType,
  text: PropTypes.string.isRequired,
  options: PropTypes.arrayOf(PropTypes.string).isRequired,
  values: PropTypes.arrayOf(PropTypes.string),
  onSelect: PropTypes.func,
};

DropdownButton.defaultProps = {
  button: Button,
  values: null,
  onSelect: () => {},
};

export default DropdownButton;
