import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import { updateQueryStringParam } from '../../utils/updateQueryParam';
import { TABS } from '../../constants/defs.json';

import DropdownButton from '../Buttons/DropdownButton';

import {
  TabBarContainer,
  TabBarList,
  TabItem,
  DropdownTab,
} from './TabBar.styles';

const TabBar = observer(({ activeTab, onTabChange }) => {
  const { configStore } = useStores();

  const changeTab = (tab, e) => {
    if (e !== undefined) {
      e.preventDefault();
    }
    onTabChange(tab);
    updateQueryStringParam('tab', tab);
  };

  const onMiscTabSelect = (tab) => {
    changeTab(tab);
  };

  const tabs = [
    <TabItem key={TABS.TAB_EXAMPLE} active={activeTab === TABS.TAB_EXAMPLE}>
      <a
        href="#"
        className="tab-link"
        onClick={changeTab.bind(this, TABS.TAB_EXAMPLE)}
      >
        <span>Home</span>
      </a>
    </TabItem>,
    <TabItem key={TABS.TAB_GROUP} active={activeTab === TABS.TAB_GROUP}>
      <a
        href="#"
        className="tab-link"
        onClick={changeTab.bind(this, TABS.TAB_GROUP)}
      >
        <span>Compare {configStore.getGroupLabel()}s</span>
      </a>
    </TabItem>,
    <TabItem key={TABS.TAB_LOCATION} active={activeTab === TABS.TAB_LOCATION}>
      <a
        href="#"
        className="tab-link"
        onClick={changeTab.bind(this, TABS.TAB_LOCATION)}
      >
        <span>Compare Locations</span>
      </a>
    </TabItem>,
    <TabItem
      key={TABS.TAB_GLOBAL_SEQUENCES}
      active={activeTab === TABS.TAB_GLOBAL_SEQUENCES}
    >
      <a
        href="#"
        className="tab-link"
        onClick={changeTab.bind(this, TABS.TAB_GLOBAL_SEQUENCES)}
      >
        <span>Global Sequencing Coverage</span>
      </a>
    </TabItem>,
    <TabItem key={TABS.TAB_ABOUT} active={activeTab === TABS.TAB_ABOUT}>
      <a
        href="#"
        className="tab-link"
        onClick={changeTab.bind(this, TABS.TAB_ABOUT)}
      >
        <span>Acknowledgements</span>
      </a>
    </TabItem>,
    <DropdownButton
      key="tab-dropdown"
      button={DropdownTab}
      text="More..."
      options={['Methods', 'Related Projects']}
      values={[TABS.TAB_METHODOLOGY, TABS.TAB_RELATED]}
      onSelect={onMiscTabSelect}
      active={[TABS.TAB_METHODOLOGY, TABS.TAB_RELATED].includes(activeTab)}
    />,
  ];

  return (
    <TabBarContainer height={tabs.length * 30}>
      <TabBarList>{tabs}</TabBarList>
    </TabBarContainer>
  );
});
TabBar.propTypes = {
  activeTab: PropTypes.string.isRequired,
  onTabChange: PropTypes.func.isRequired,
};
TabBar.defaultProps = {};

export default TabBar;
