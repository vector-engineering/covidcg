import React, { useState } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';

const LegendAndButtonWrapper = styled.div`
  display: flex;
  flex-direction: row;
  align-items: flex-start;
  margin-top: 12px;
  margin-bottom: 6px;
  margin-left: 24px;
  margin-right: 24px;
`;

const LegendContainer = styled.div`
  width: 100%;
  padding: 4px 6px;
  transition: max-height 0.3s cubic-bezier(0.22, 1, 0.36, 1);
  align-self: center;
  border-radius: 2px;
  max-height: ${({ maxHeight, collapsed }) => {
    if (collapsed) {
      return '1px';
    }
    return maxHeight ? `${maxHeight}` : '100px';
  }};
  overflow-y: scroll;
`;

const CollapseButton = styled.button`
  border: none;
  cursor: pointer;
  color: #ccc;
  transition: color 0.2s cubic-bezier(0.22, 1, 0.36, 1);
  background-color: rgba(0, 0, 0, 0);
  outline: none;
  width: 32px;
  font-size: 20px;

  &:hover {
    color: black;
  }
`;

const Title = styled.div`
  padding-left: 12px;
  display: flex;
`;

const AccordionWrapper = ({ children, maxHeight, title, defaultCollapsed }) => {
  const [collapsed, setCollapsed] = useState(defaultCollapsed);
  return (
    <div style={{ width: '100%' }}>
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
      <LegendAndButtonWrapper>
        <LegendContainer maxHeight={maxHeight} collapsed={collapsed}>
          {children}
        </LegendContainer>
      </LegendAndButtonWrapper>
    </div>
  );
};

AccordionWrapper.propTypes = {
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node,
  ]),
  maxHeight: PropTypes.number.isRequired,
  title: PropTypes.string.isRequired,
  defaultCollapsed: PropTypes.bool.isRequired,
};

export default AccordionWrapper;
