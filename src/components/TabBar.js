import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
// import { useStores } from '../stores/connect';
import styled from 'styled-components';
// import _ from 'underscore';

// import { intToISO } from '../utils/date';

const TabBarContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: stretch;

  padding-top: 5px;
  padding-left: 10px;
  padding-right: 10px;

  width: 100%;
  min-height: 30px;
  margin-bottom: 10px;
  background-color: #eee;
`;

const TabBarList = styled.div`
  display: flex;
  flex-direction: row;
  align-items: stretch;
  justify-content: flex-start;
  flex-grow: 1;
`;

const TabItem = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  width: 180px;

  a.tab-link {
    display: flex;
    flex-direction: column;
    align-items: stretch;
    justify-content: center;

    flex-grow: 1;
    text-decoration: none;
    color: #000;
    text-align: center;
    border-radius: 10px 10px 0px 0px;
    background-color: ${({ active }) => (active ? '#fff' : 'transparent')};
    transition: 0.1s all ease-in-out;

    &:hover {
      color: #666;
      background-color: ${({ active }) => (active ? '#fff' : '#f8f8f8')};
    }
    &:active {
      color: #888;
    }

    span {
      border-right: ${({ active }) => (active ? 'none' : '1px solid #ccc')};
    }
  }
`;
TabItem.defaultProps = {
  active: false,
};

// const StatusGroup = styled.div`
//   display: flex;
//   flex-direction: row;
//   align-items: stretch;
//   justify-content: flex-start;
//   flex-grow: 1;

//   border-bottom: 1px solid #aaa;
//   // box-shadow: 0px 0px 3px #ccc;
// `;

// const RowStatus = styled.div`
//   display: flex;
//   flex-direction: column;
//   align-items: flex-start;
//   justify-content: flex-start;

//   flex-grow: ${({ grow }) => grow};

//   padding: 4px 10px;
//   line-height: normal;
//   border-right: 1px solid #ccc;

//   .row-status-header {
//     .title {
//       font-size: 0.9em;
//       margin-right: 5px;
//     }
//     .subtitle {
//       font-size: 0.8em;
//       font-weight: normal;
//     }
//   }

//   .content {
//     font-weight: normal;
//     font-size: 0.8em;
//     text-overflow: ellipsis;
//     white-space: nowrap;
//     width: 100%;
//     overflow: hidden;
//   }
// `;
// RowStatus.defaultProps = {
//   width: 200,
// };

const TabBar = observer(({ activeTab, onTabChange }) => {
  // const { dataStore } = useStores();

  const changeTab = (tab, e) => {
    e.preventDefault();
    onTabChange(tab);
  };

  // let numSequencesBeforeMetadataFilteringText = '';
  // if (
  //   dataStore.selectedRows.length !==
  //   dataStore.numSequencesBeforeMetadataFiltering
  // ) {
  //   numSequencesBeforeMetadataFilteringText = (
  //     <span className="content">
  //       ({dataStore.numSequencesBeforeMetadataFiltering} before metadata
  //       filtering)
  //     </span>
  //   );
  // }

  // let groupName;
  // if (dataStore.groupKey === 'lineage') {
  //   groupName = 'Lineage';
  // } else if (dataStore.groupKey === 'clade') {
  //   groupName = 'Clade';
  // } else if (dataStore.groupKey === 'snp') {
  //   if (dataStore.dnaOrAa === 'dna') {
  //     groupName = 'NT SNP';
  //   } else {
  //     groupName = 'AA SNP';
  //   }
  // }

  // let selectedGroups = '';
  // if (dataStore.selectedGroups.length > 0) {
  //   selectedGroups = _.pluck(dataStore.selectedGroups, 'group').join(', ');
  //   selectedGroups = <span className="content">({selectedGroups})</span>;
  // }

  // let selectedDates = 'No date range selected';
  // if (dataStore.dateRange[0] !== -1 && dataStore.dateRange[1] !== -1) {
  //   selectedDates =
  //     intToISO(dataStore.dateRange[0]) +
  //     ' – ' +
  //     intToISO(dataStore.dateRange[1] - 86400000);
  // }

  return (
    <TabBarContainer>
      <TabBarList>
        <TabItem active={activeTab === 'group'}>
          <a
            href="#"
            className="tab-link"
            onClick={changeTab.bind(this, 'group')}
          >
            <span>Main</span>
          </a>
        </TabItem>
        <TabItem active={activeTab === 'location'}>
          <a
            href="#"
            className="tab-link"
            onClick={changeTab.bind(this, 'location')}
          >
            <span>Compare locations</span>
          </a>
        </TabItem>
        <TabItem active={activeTab === 'about'}>
          <a
            href="#"
            className="tab-link"
            onClick={changeTab.bind(this, 'about')}
          >
            <span>About</span>
          </a>
        </TabItem>
      </TabBarList>
      {/* <StatusGroup>
        <RowStatus grow={1}>
          <div className="row-status-header">
            <span className="title">Sequences</span>
            <span className="subtitle">
              ({dataStore.selectedRows.length} selected)
            </span>
          </div>
          {numSequencesBeforeMetadataFilteringText}
        </RowStatus>
        <RowStatus grow={2}>
          <div className="row-status-header">
            <span className="title">{groupName} Selection</span>
            <span className="subtitle">
              ({dataStore.selectedGroups.length} selected)
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
      </StatusGroup> */}
    </TabBarContainer>
  );
});
TabBar.propTypes = {
  activeTab: PropTypes.string.isRequired,
  onTabChange: PropTypes.func.isRequired,
};
TabBar.defaultProps = {};

export default TabBar;
