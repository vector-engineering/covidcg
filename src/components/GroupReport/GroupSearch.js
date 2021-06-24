import React, { useMemo, useState, useEffect } from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import SelectedGroups from './SelectedGroups';

import {
  GroupSearchContainer,
  GroupSearchHeader,
  GroupSearchTitle,
  StyledDropdownTreeSelect,
} from './GroupSearch.styles';

const GroupSearch = observer(() => {
  const { groupDataStore } = useStores();
  const [state, setState] = useState({
    data: [],
  });

  useEffect(() => {
    setState({
      ...state,
      data: groupDataStore.groupSelectTree[groupDataStore.activeGroupType],
    });
  }, [groupDataStore.groupSelectTree]);

  const treeSelectOnChange = (_currentNode, selectedNodes) => {
    // console.log('onChange::', _currentNode, selectedNodes);
    groupDataStore.updateSelectedGroups(
      selectedNodes.map((node) => {
        return node.value;
      })
    );
  };
  // const treeSelectOnAction = (node, action) => {
  //   console.log('onAction::', action, node);
  // };

  const dropdownContainer = useMemo(
    () => (
      <StyledDropdownTreeSelect
        mode={'multiSelect'}
        data={state.data}
        className="geo-dropdown-tree-select"
        clearSearchOnChange={false}
        keepTreeOnSearch={true}
        keepChildrenOnSearch={true}
        showPartiallySelected={false}
        inlineSearchInput={false}
        texts={{
          placeholder: 'Search...',
          noMatches: 'No matches found',
        }}
        onChange={treeSelectOnChange}
        // onAction={treeSelectOnAction}
        // onNodeToggle={treeSelectOnNodeToggleCurrentNode}
      />
    ),
    [state.data]
  );

  return (
    <GroupSearchContainer>
      <SelectedGroups />
      <GroupSearchHeader>
        <GroupSearchTitle>
          Select {groupDataStore.getActiveGroupTypePrettyName()}s
        </GroupSearchTitle>
      </GroupSearchHeader>
      {dropdownContainer}
    </GroupSearchContainer>
  );
});

export default GroupSearch;
