import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import Button from '../Buttons/Button';

import {
  SelectedGroupsContainer,
  SelectedGroupsList,
  SelectedGroupItemContainer,
  SelectedGroupItemTitle,
  NoGroupsSelectedContainer,
} from './SelectedGroups.styles';

const SelectedGroupItem = ({ name, onClick }) => {
  return (
    <SelectedGroupItemContainer>
      <SelectedGroupItemTitle>{name}</SelectedGroupItemTitle>
      <Button onClick={onClick}>x</Button>
    </SelectedGroupItemContainer>
  );
};
SelectedGroupItem.propTypes = {
  name: PropTypes.string,
  onClick: PropTypes.func,
};

const SelectedGroups = observer(() => {
  const { groupDataStore } = useStores();

  const removeGroup = (groupNameToRemove) => {
    groupDataStore.updateSelectedGroups(
      groupDataStore.selectedGroups.slice().filter((groupName) => {
        return groupName != groupNameToRemove;
      })
    );
  };

  const selectedGroupItems = [];
  groupDataStore.selectedGroups.forEach((group) => {
    selectedGroupItems.push(
      <SelectedGroupItem
        key={`report-selected-group-${group}`}
        name={group}
        onClick={removeGroup.bind(this, group)}
      />
    );
  });

  const selectedGroupsList = (
    <SelectedGroupsList>{selectedGroupItems}</SelectedGroupsList>
  );

  const noGroupsSelected = (
    <NoGroupsSelectedContainer>
      No {groupDataStore.getActiveGroupTypePrettyName()}s selected
    </NoGroupsSelectedContainer>
  );

  return (
    <SelectedGroupsContainer>
      {/* Display separate element if no groups are selected */}
      {groupDataStore.selectedGroups.length == 0 && noGroupsSelected}
      {groupDataStore.selectedGroups.length > 0 && selectedGroupsList}
    </SelectedGroupsContainer>
  );
});

export default SelectedGroups;
