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

const DropdownButton = ({ text, options, onSelect }) => {
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
      handleToggle();
    }
  };
  // const handleBlur = () => {};

  const handleOption = (e) => {
    onSelect(e.target.dataset.option);
    handleToggle();
  };

  // Build option anchor links
  const optionElements = [];
  options.forEach((option) => {
    optionElements.push(
      <DropdownItem
        key={option}
        className="dropdown-item"
        href="#"
        onClick={handleOption}
        data-option={option}
      >
        {option}
      </DropdownItem>
    );
  });

  return (
    <DropdownContainer>
      <Button
        id={state.id}
        data-toggle="dropdown"
        aria-haspopup="true"
        aria-expanded={state.expanded}
        onClick={handleToggle}
        // onFocus={handleFocus}
        onBlur={handleBlur}
        sticky={true}
      >
        {text}
        <span className="caret"></span>
      </Button>
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
  text: PropTypes.string.isRequired,
  options: PropTypes.arrayOf(PropTypes.string).isRequired,
  onSelect: PropTypes.func,
};

DropdownButton.defaultProps = {
  onSelect: () => {},
};

export default DropdownButton;
