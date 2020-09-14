import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';

import { updateQueryStringParam } from '../../utils/updateQueryParam';

const TabIndicatorBG = styled.a`
  display: inline-flex;
  flex-direction: column;
  align-items: stretch;

  width: 130px;

  background-color: #ddd;
  padding: 2px 3px 0px 3px;
  border: 1px solid #ccc;
  border-bottom: none;
  margin-top: -1px;

  text-decoration: none;
  color: #444;

  &:hover {
    span {
      background-color: #f8f8f8;
    }
  }
`;

const TabIndicatorFG = styled.span`
  font-size: 0.85em;
  font-weight: 500;
  background-color: #fff;
  border-radius: 5px 5px 0px 0px;
  text-align: center;
  transition: 0.1s all ease-in-out;
`;

const TabIndicatorWrapper = ({ tab, children }) => {
  const { UIStore } = useStores();

  const handleClick = (e) => {
    e.preventDefault();

    if (tab !== null) {
      UIStore.setActiveTab(tab);
      updateQueryStringParam('tab', tab);
    }
  };

  return (
    <TabIndicatorBG disabled={tab === null} href="#" onClick={handleClick}>
      <TabIndicatorFG>{children}</TabIndicatorFG>
    </TabIndicatorBG>
  );
};
TabIndicatorWrapper.propTypes = {
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node,
  ]).isRequired,
  tab: PropTypes.string,
};
TabIndicatorWrapper.defaultProps = {
  tab: null,
};

export default TabIndicatorWrapper;
