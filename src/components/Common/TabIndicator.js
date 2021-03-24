import React from 'react';
import PropTypes from 'prop-types';
import { useStores } from '../../stores/connect';

import { updateQueryStringParam } from '../../utils/updateQueryParam';

import { TabIndicatorBG, TabIndicatorFG } from './TabIndicator.styles';

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
