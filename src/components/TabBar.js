import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../stores/connect';
import styled from 'styled-components';
// import _ from 'underscore';

// import { intToISO } from '../utils/date';
// import { GROUP_KEYS, DNA_OR_AA } from '../constants/config';
import { TABS } from '../constants/UI';
import { updateQueryStringParam } from '../utils/updateQueryParam';

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

  max-width: 160px;

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
  const { configStore } = useStores();

  const changeTab = (tab, e) => {
    e.preventDefault();
    onTabChange(tab);
    updateQueryStringParam('tab', tab);
  };

  // let numSequencesBeforeMetadataFilteringText = '';
  // if (
  //   dataStore.selectedAccessionIds.length !==
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
  // if (configStore.groupKey === GROUP_KEYS.GROUP_LINEAGE) {
  //   groupName = 'Lineage';
  // } else if (configStore.groupKey === GROUP_KEYS.GROUP_CLADE) {
  //   groupName = 'Clade';
  // } else if (configStore.groupKey === GROUP_KEYS.GROUP_SNV) {
  //   if (configStore.dnaOrAa === DNA_OR_AA.DNA) {
  //     groupName = 'NT SNP';
  //   } else {
  //     groupName = 'AA SNP';
  //   }
  // }

  // let selectedGroups = '';
  // if (configStore.selectedGroups.length > 0) {
  //   selectedGroups = _.pluck(configStore.selectedGroups, 'group').join(', ');
  //   selectedGroups = <span className="content">({selectedGroups})</span>;
  // }

  // let selectedDates = 'No date range selected';
  // if (configStore.dateRange[0] !== -1 && configStore.dateRange[1] !== -1) {
  //   selectedDates =
  //     intToISO(configStore.dateRange[0]) +
  //     ' – ' +
  //     intToISO(configStore.dateRange[1] - 86400000);
  // }

  return (
    <TabBarContainer>
      <TabBarList>
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
        <TabItem active={activeTab === TABS.TAB_EXAMPLE}>
          <a
            href="#"
            className="tab-link"
            onClick={changeTab.bind(this, TABS.TAB_EXAMPLE)}
          >
            <span>Analyses</span>
          </a>
        </TabItem>
        <TabItem active={activeTab === TABS.TAB_GLOBAL_SEQUENCES}>
          <a
            href="#"
            className="tab-link"
            onClick={changeTab.bind(this, TABS.TAB_GLOBAL_SEQUENCES)}
          >
            <span>Sequencing Efforts</span>
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
        <TabItem active={activeTab === TABS.TAB_METHODOLOGY}>
          <a
            href="#"
            className="tab-link"
            onClick={changeTab.bind(this, TABS.TAB_METHODOLOGY)}
          >
            <span>Methods</span>
          </a>
        </TabItem>
        <TabItem active={activeTab === TABS.TAB_RELATED}>
          <a
            href="#"
            className="tab-link"
            onClick={changeTab.bind(this, TABS.TAB_RELATED)}
          >
            <span>Related Projects</span>
          </a>
        </TabItem>
      </TabBarList>
      {/* <StatusGroup>
        <RowStatus grow={1}>
          <div className="row-status-header">
            <span className="title">Sequences</span>
            <span className="subtitle">
              ({dataStore.selectedAccessionIds.length} selected)
            </span>
          </div>
          {numSequencesBeforeMetadataFilteringText}
        </RowStatus>
        <RowStatus grow={2}>
          <div className="row-status-header">
            <span className="title">{groupName} Selection</span>
            <span className="subtitle">
              ({configStore.selectedGroups.length} selected)
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
