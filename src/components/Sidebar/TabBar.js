import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import { TABS } from '../../constants/defs.json';

import DropdownButton from '../Buttons/DropdownButton';

import {
  TabBarContainer,
  TabBarList,
  TabItem,
  DropdownTab,
} from './TabBar.styles';

import { config } from '../../config';

const TabBar = observer(({ activeTab, onTabChange }) => {
  const { configStore } = useStores();

  const changeTab = (tab, e) => {
    if (e !== undefined) {
      e.preventDefault();
    }
    onTabChange(tab);
  };

  const onMiscTabSelect = (tab) => {
    changeTab(tab);
  };

  // The base tabs are the Example, Compare Groups, and Compare Locations tabs
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
  ];

  // Lineage Reports tab should go second when it functions
  if (config.virus === 'sars2') {
    tabs.push(
      <TabItem
        key={TABS.TAB_GROUP_REPORT}
        active={activeTab === TABS.TAB_GROUP_REPORT}
      >
        <a
          href="#"
          className="tab-link"
          onClick={changeTab.bind(this, TABS.TAB_GROUP_REPORT)}
        >
          <span>Lineage Reports</span>
        </a>
      </TabItem>
    );
  }

  // Add Compare Groups/Locations
  tabs.push(
    <TabItem
      key={TABS.TAB_COMPARE_GROUPS}
      active={activeTab === TABS.TAB_COMPARE_GROUPS}
    >
      <a
        href="#"
        className="tab-link"
        onClick={changeTab.bind(this, TABS.TAB_COMPARE_GROUPS)}
      >
        <span>Compare {configStore.getGroupLabel()}s</span>
      </a>
    </TabItem>,
    <TabItem
      key={TABS.TAB_COMPARE_LOCATIONS}
      active={activeTab === TABS.TAB_COMPARE_LOCATIONS}
    >
      <a
        href="#"
        className="tab-link"
        onClick={changeTab.bind(this, TABS.TAB_COMPARE_LOCATIONS)}
      >
        <span>Compare Locations</span>
      </a>
    </TabItem>
  );

  // Add virus specific tabs
  if (config.virus === 'sars2') {
    tabs.push(
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
      </TabItem>
    );
  }

  // Add About Tab and DropdownButton last
  tabs.push(
    <TabItem key={TABS.TAB_ABOUT} active={activeTab === TABS.TAB_ABOUT}>
      <a
        href="#"
        className="tab-link"
        onClick={changeTab.bind(this, TABS.TAB_ABOUT)}
      >
        {config.virus === 'sars2' && <span>About COVID CG</span>}
        {config.virus === 'rsv' && <span>About RSV CG</span>}
      </a>
    </TabItem>
  );

  tabs.push(
    <DropdownButton
      key="tab-dropdown"
      button={DropdownTab}
      text="More..."
      options={['Methods', 'Related Projects']}
      values={[TABS.TAB_METHODOLOGY, TABS.TAB_RELATED]}
      onSelect={onMiscTabSelect}
      active={[TABS.TAB_METHODOLOGY, TABS.TAB_RELATED].includes(activeTab)}
    />
  );

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
