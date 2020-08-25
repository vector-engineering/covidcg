import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';

const TabIndicatorBG = styled.span`
  display: inline-flex;
  flex-direction: column;
  align-items: stretch;

  width: 130px;

  background-color: #ddd;
  padding: 2px 3px 0px 3px;
  border: 1px solid #ccc;
  border-bottom: none;
  margin-top: -1px;
`;

const TabIndicatorFG = styled.span`
  font-size: 0.85em;
  font-weight: 500;
  background-color: #fff;
  border-radius: 5px 5px 0px 0px;
  text-align: center;
`;

const TabIndicatorWrapper = ({ children }) => {
  return (
    <TabIndicatorBG>
      <TabIndicatorFG>{children}</TabIndicatorFG>
    </TabIndicatorBG>
  );
};
TabIndicatorWrapper.propTypes = {
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node,
  ]).isRequired,
};

export default TabIndicatorWrapper;
