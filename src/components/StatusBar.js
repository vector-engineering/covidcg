import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../stores/connect';
import styled from 'styled-components';
import _ from 'underscore';

import { intToISO } from '../utils/date';

const StatusBarContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: stretch;

  width: 100%;
  min-height: 40px;
  margin-bottom: 10px;
  background-color: #f8f8f8;
`;

const TabBar = styled.div`
  display: flex;
  flex-direction: row;
  align-items: stretch;
  justify-content: flex-start;
`;

const TabItem = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  width: 150px;
  background-color: ${({ active }) => (active ? '#fff' : 'transparent')};
  ${({ active }) => {
    if (active) {
      return `
        border-width: 1px 1px 0px 1px;
        border-color: #aaa;
        border-style: solid;
      `;
    } else {
      return `
        border-right: 1px solid #ddd;
        border-bottom: 1px solid #aaa;
      `;
    }
  }}
  // border-radius: 5px 5px 0px 0px;

  &:hover {
    background-color: ${({ active }) => (active ? '#fff' : '#eee')};
  }

  a.tab-link {
    padding: 10px 8px 0px 8px;
    flex-grow: 1;
    text-decoration: none;
    color: #000;
    text-align: center;

    &:hover {
      color: #666;
    }
    &:active {
      color: #888;
    }
  }
`;
TabItem.defaultProps = {
  active: false,
};

const StatusGroup = styled.div`
  display: flex;
  flex-direction: row;
  align-items: stretch;
  justify-content: flex-start;
  flex-grow: 1;

  border-bottom: 1px solid #aaa;
  // box-shadow: 0px 0px 3px #ccc;
`;

const RowStatus = styled.div`
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  justify-content: flex-start;

  flex-grow: ${({ grow }) => grow};

  padding: 4px 10px;
  line-height: normal;
  border-right: 1px solid #ccc;

  .row-status-header {
    .title {
      font-size: 0.9em;
      margin-right: 5px;
    }
    .subtitle {
      font-size: 0.8em;
      font-weight: normal;
    }
  }

  .content {
    font-weight: normal;
    font-size: 0.8em;
    text-overflow: ellipsis;
    white-space: nowrap;
    width: 100%;
    overflow: hidden;
  }
`;
RowStatus.defaultProps = {
  width: 200,
};

const StatusBar = observer(({ activeTab, onTabChange }) => {
  const { covidStore } = useStores();

  const changeTab = (tab, e) => {
    e.preventDefault();
    onTabChange(tab);
  };

  let numSequencesBeforeMetadataFilteringText = '';
  if (
    covidStore.selectedRows.length !==
    covidStore.numSequencesBeforeMetadataFiltering
  ) {
    numSequencesBeforeMetadataFilteringText = (
      <span className="content">
        ({covidStore.numSequencesBeforeMetadataFiltering} before metadata
        filtering)
      </span>
    );
  }

  let groupName;
  if (covidStore.groupKey === 'lineage') {
    groupName = 'Lineage';
  } else if (covidStore.groupKey === 'clade') {
    groupName = 'Clade';
  } else if (covidStore.groupKey === 'snp') {
    if (covidStore.dnaOrAa === 'dna') {
      groupName = 'NT SNP';
    } else {
      groupName = 'AA SNP';
    }
  }

  let selectedGroups = '';
  if (covidStore.selectedGroups.length > 0) {
    selectedGroups = _.pluck(covidStore.selectedGroups, 'group').join(', ');
    selectedGroups = <span className="content">({selectedGroups})</span>;
  }

  let selectedDates = 'No date range selected';
  if (covidStore.dateRange[0] !== -1 && covidStore.dateRange[1] !== -1) {
    selectedDates =
      intToISO(covidStore.dateRange[0]) +
      ' – ' +
      intToISO(covidStore.dateRange[1] - 86400000);
  }

  return (
    <StatusBarContainer>
      <TabBar>
        <TabItem active={activeTab === 'group'}>
          <a
            href="#"
            className="tab-link"
            onClick={changeTab.bind(this, 'group')}
          >
            Main
          </a>
        </TabItem>
        <TabItem active={activeTab === 'location'}>
          <a
            href="#"
            className="tab-link"
            onClick={changeTab.bind(this, 'location')}
          >
            Compare locations
          </a>
        </TabItem>
        <TabItem active={activeTab === 'about'}>
          <a
            href="#"
            className="tab-link"
            onClick={changeTab.bind(this, 'about')}
          >
            About
          </a>
        </TabItem>
      </TabBar>
      <StatusGroup>
        <RowStatus grow={1}>
          <div className="row-status-header">
            <span className="title">Sequences</span>
            <span className="subtitle">
              ({covidStore.selectedRows.length} selected)
            </span>
          </div>
          {numSequencesBeforeMetadataFilteringText}
        </RowStatus>
        <RowStatus grow={2}>
          <div className="row-status-header">
            <span className="title">{groupName} Selection</span>
            <span className="subtitle">
              ({covidStore.selectedGroups.length} selected)
            </span>
          </div>
          {selectedGroups}
        </RowStatus>
        <RowStatus grow={1}>
          <div className="row-status-header">
            <span className="title">Date Selection</span>
          </div>
          <span className="content">{selectedDates}</span>
        </RowStatus>
      </StatusGroup>
    </StatusBarContainer>
  );
});
StatusBar.propTypes = {
  activeTab: PropTypes.string.isRequired,
  onTabChange: PropTypes.func.isRequired,
};
StatusBar.defaultProps = {};

export default StatusBar;
