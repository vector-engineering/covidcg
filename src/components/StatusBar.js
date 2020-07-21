import React from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../stores/connect';
import styled from 'styled-components';
import _ from 'underscore';

const StatusBarContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: stretch;

  width: 100%;
  min-height: 40px;
  margin-bottom: 10px;
  border-bottom: 1px solid #aaa;
  // box-shadow: 0px 0px 3px #ccc;

  background-color: #f8f8f8;
`;

const RowStatus = styled.div`
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  justify-content: flex-start;

  width: ${({ width }) => width}px;

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

const StatusBar = observer(() => {
  const { covidStore } = useStores();

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

  return (
    <StatusBarContainer>
      <RowStatus width={200}>
        <div className="row-status-header">
          <span className="title">Sequences</span>
          <span className="subtitle">
            ({covidStore.selectedRows.length} selected)
          </span>
        </div>
        {numSequencesBeforeMetadataFilteringText}
      </RowStatus>
      <RowStatus width={300}>
        <div className="row-status-header">
          <span className="title">{groupName} Selection</span>
          <span className="subtitle">
            ({covidStore.selectedGroups.length} selected)
          </span>
        </div>
        {selectedGroups}
      </RowStatus>
    </StatusBarContainer>
  );
});

export default StatusBar;
