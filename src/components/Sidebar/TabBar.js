import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import { config } from '../../config';
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
  ];

  if (config['show_reports_tab']) {
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

  if (config['show_global_sequencing_tab']) {
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

  tabs.push(
    <TabItem key={TABS.TAB_ABOUT} active={activeTab === TABS.TAB_ABOUT}>
      <a
        href="#"
        className="tab-link"
        onClick={changeTab.bind(this, TABS.TAB_ABOUT)}
      >
        {config.virus === 'sars2' && <span>About COVID CG</span>}
        {config.virus === 'rsv' && <span>About RSV PathMut</span>}
      </a>
    </TabItem>
  );

  const misc_tabs = [];
  const misc_tab_titles = [];
  if (config['show_methods_tab']) {
    misc_tabs.push(TABS.TAB_METHODOLOGY);
    misc_tab_titles.push('Methods');
  }
  if (config['show_related_projects_tab']) {
    misc_tabs.push(TABS.TAB_RELATED);
    misc_tab_titles.push('Related Projects');
  }

  if (misc_tabs.length > 0) {
    tabs.push(
      <DropdownButton
        key="tab-dropdown"
        button={DropdownTab}
        text="More..."
        options={misc_tab_titles}
        values={misc_tabs}
        onSelect={onMiscTabSelect}
        active={misc_tabs.includes(activeTab)}
      />
    );
  }

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
