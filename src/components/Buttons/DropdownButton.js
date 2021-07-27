import React, { useState } from 'react';
import PropTypes from 'prop-types';

import Button from './Button';

import {
  DropdownContainer,
  DropdownMenu,
  DropdownItem,
} from './DropdownButton.styles';

import { uniqueId } from '../../utils/func';

const DropdownButton = ({
  button,
  text,
  options,
  values,
  onSelect,
  direction,
  ...other
}) => {
  const [state, setState] = useState({
    expanded: false,
    id: uniqueId('dropdown_'),
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
        direction={direction}
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
  direction: PropTypes.string,
};

DropdownButton.defaultProps = {
  button: Button,
  values: null,
  onSelect: () => {},
  direction: 'right',
};

export default DropdownButton;
