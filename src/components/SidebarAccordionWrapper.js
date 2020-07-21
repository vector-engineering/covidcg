import React, { useState } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';

const AccordionContainer = styled.div`
  width: 100%;
  margin-bottom: 3px;
  border-top: 1px solid #aaa;
  padding-top: 5px;
`;

const CollapseContent = styled.div`
  width: 100%;
  transition: max-height 0.15s cubic-bezier(0.22, 1, 0.36, 1);
  align-self: center;
  border-radius: 2px;
  max-height: ${({ collapsed }) => {
    if (collapsed) {
      return '1px';
    } else {
      return '100%';
    }
  }};
  overflow-y: scroll;
`;

const CollapseButton = styled.button`
  border: none;
  cursor: pointer;
  color: #aaa;
  transition: color 0.15s cubic-bezier(0.22, 1, 0.36, 1);
  background-color: rgba(0, 0, 0, 0);
  outline: none;
  width: 16px;
  font-size: 16px;
  font-family: monospace;
  padding: 0px;
  margin-left: 10px;
  margin-right: 3px;

  &:hover {
    color: black;
  }
`;

const Title = styled.div`
  display: flex;
`;

const SidebarAccordionWrapper = ({
  children,
  maxHeight,
  title,
  defaultCollapsed,
}) => {
  const [collapsed, setCollapsed] = useState(defaultCollapsed);
  return (
    <AccordionContainer>
      <Title>
        <CollapseButton
          onClick={() => {
            setCollapsed(!collapsed);
          }}
        >
          {collapsed ? '+' : '-'}
        </CollapseButton>
        {title}
      </Title>
      <CollapseContent maxHeight={maxHeight} collapsed={collapsed}>
        {children}
      </CollapseContent>
    </AccordionContainer>
  );
};

SidebarAccordionWrapper.propTypes = {
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node,
  ]),
  maxHeight: PropTypes.string.isRequired,
  title: PropTypes.string.isRequired,
  defaultCollapsed: PropTypes.bool.isRequired,
};

export default SidebarAccordionWrapper;
