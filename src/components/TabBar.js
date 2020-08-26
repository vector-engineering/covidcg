import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../stores/connect';
import styled from 'styled-components';
import { updateQueryStringParam } from '../utils/updateQueryParam';

import DropdownButton from './Buttons/DropdownButton';

import { TABS } from '../constants/UI';

const TabBarContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: stretch;

  padding-top: 5px;

  width: 100%;
  min-height: 30px;
  background-color: #eee;
`;

const TabBarList = styled.div`
  display: flex;
  flex-direction: row;
  align-items: stretch;
  justify-content: flex-start;
  flex-grow: 1;

  padding-left: 10px;
  padding-right: 10px;
`;

const TabItem = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  max-width: 250px;

  a.tab-link {
    display: flex;
    flex-direction: column;
    align-items: stretch;
    justify-content: center;

    flex-grow: 1;
    text-decoration: none;
    color: ${({ active }) => (active ? '#000' : '#888')};
    text-align: center;
    border-radius: 10px 10px 0px 0px;
    background-color: ${({ active }) => (active ? '#fff' : 'transparent')};
    transition: 0.1s all ease-in-out;

    padding: 0px 15px;

    &:hover {
      color: #666;
      background-color: ${({ active }) => (active ? '#fff' : '#f8f8f8')};
    }
    &:active {
      color: #888;
    }

    span {
      // border-right: ${({ active }) => (active ? 'none' : '1px solid #ccc')};
    }
  }
`;
TabItem.defaultProps = {
  active: false,
};

const DropdownTab = styled.button`
  display: flex;
  flex-direction: row;
  align-items: center;

  height: 30px;
  padding: 0px 15px;
  padding-right: 20px;

  border: none;
  outline: none;

  background-color: ${({ active }) => (active ? '#fff' : 'transparent')};
  line-height: normal;
  text-align: center;
  font: 14px 'Helvetica Neue', Helvetica, Arial, sans-serif;
  font-weight: 500;
  color: ${({ active }) => (active ? '#000' : '#888')};

  border-radius: 10px 10px 0px 0px;
  transition: 0.1s all ease-in-out;

  .caret {
    height: 5px;

    &:after {
      left: 1px;
      top: 0px;
      border-top: 5px solid #888;
    }
  }
`;
DropdownTab.defaultProps = {
  active: false,
};

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

  return (
    <TabBarContainer>
      <TabBarList>
        <TabItem active={activeTab === TABS.TAB_EXAMPLE}>
          <a
            href="#"
            className="tab-link"
            onClick={changeTab.bind(this, TABS.TAB_EXAMPLE)}
          >
            <span>Example Analyses</span>
          </a>
        </TabItem>
        <TabItem active={activeTab === TABS.TAB_GROUP}>
          <a
            href="#"
            className="tab-link"
            onClick={changeTab.bind(this, TABS.TAB_GROUP)}
          >
            <span>Compare {configStore.getGroupLabel()}s</span>
          </a>
        </TabItem>
        <TabItem active={activeTab === TABS.TAB_LOCATION}>
          <a
            href="#"
            className="tab-link"
            onClick={changeTab.bind(this, TABS.TAB_LOCATION)}
          >
            <span>Compare Locations</span>
          </a>
        </TabItem>
        <TabItem active={activeTab === TABS.TAB_GLOBAL_SEQUENCES}>
          <a
            href="#"
            className="tab-link"
            onClick={changeTab.bind(this, TABS.TAB_GLOBAL_SEQUENCES)}
          >
            <span>Global Sequencing Coverage</span>
          </a>
        </TabItem>
        <TabItem active={activeTab === TABS.TAB_ABOUT}>
          <a
            href="#"
            className="tab-link"
            onClick={changeTab.bind(this, TABS.TAB_ABOUT)}
          >
            <span>Acknowledgements</span>
          </a>
        </TabItem>
        <DropdownButton
          button={DropdownTab}
          text="More..."
          options={['Methods', 'Related Projects']}
          values={[TABS.TAB_METHODOLOGY, TABS.TAB_RELATED]}
          onSelect={onMiscTabSelect}
          active={[TABS.TAB_METHODOLOGY, TABS.TAB_RELATED].includes(activeTab)}
        />
      </TabBarList>
    </TabBarContainer>
  );
});
TabBar.propTypes = {
  activeTab: PropTypes.string.isRequired,
  onTabChange: PropTypes.func.isRequired,
};
TabBar.defaultProps = {};

export default TabBar;
