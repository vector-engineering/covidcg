import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import {
  SelectedReportGroupTitle,
  SelectedReportGroupsContainer,
  SelectedReportGroupsList,
  SelectedReportGroupItemContainer,
  SelectedReportGroupItemTitle,
  NoReportGroupsSelectedContainer,
  SelectedReportGroupsButton,
} from './SelectedReportGroups.styles';

const SelectedGroupItem = ({ name, onClick }) => {
  return (
    <SelectedReportGroupItemContainer>
      <SelectedReportGroupsButton title="Deselect" onClick={onClick}>
        x
      </SelectedReportGroupsButton>
      <SelectedReportGroupItemTitle>{name}</SelectedReportGroupItemTitle>
    </SelectedReportGroupItemContainer>
  );
};
SelectedGroupItem.propTypes = {
  name: PropTypes.string,
  onClick: PropTypes.func,
};

const SelectedReportGroups = observer(() => {
  const { groupDataStore } = useStores();

  const removeGroup = (groupNameToRemove) => {
    groupDataStore.applyPendingChanges({
      selectedReportGroups: groupDataStore.selectedReportGroups
        .slice()
        .filter((groupName) => {
          return groupName != groupNameToRemove;
        }),
    });
  };

  const selectedGroupItems = [];
  groupDataStore.selectedReportGroups.forEach((group) => {
    selectedGroupItems.push(
      <SelectedGroupItem
        key={`report-selected-group-${group}`}
        name={group}
        onClick={removeGroup.bind(this, group)}
      />
    );
  });

  const selectedReportGroupsList = (
    <SelectedReportGroupsList>{selectedGroupItems}</SelectedReportGroupsList>
  );

  const noGroupsSelected = (
    <NoReportGroupsSelectedContainer>
      No {groupDataStore.getActiveReportGroupTypePrettyName()}s selected
    </NoReportGroupsSelectedContainer>
  );

  return (
    <SelectedReportGroupsContainer>
      <SelectedReportGroupTitle>Selected Lineages</SelectedReportGroupTitle>
      {/* Display separate element if no groups are selected */}
      {groupDataStore.selectedReportGroups.length == 0 && noGroupsSelected}
      {groupDataStore.selectedReportGroups.length > 0 &&
        selectedReportGroupsList}
    </SelectedReportGroupsContainer>
  );
});

export default SelectedReportGroups;
