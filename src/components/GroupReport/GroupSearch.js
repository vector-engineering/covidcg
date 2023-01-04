import React, { useMemo, useState, useEffect } from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import StyledDropdownTreeSelect from '../Common/StyledDropdownTreeSelect';

import {
  GroupSearchContainer,
  GroupSearchHeader,
  GroupSearchTitle,
} from './GroupSearch.styles';

const GroupSearch = observer(() => {
  const { groupDataStore, configStore } = useStores();
  const [state, setState] = useState({
    data: [],
  });

  useEffect(() => {
    let data =
      groupDataStore.reportGroupSelectTree[
        groupDataStore.activeReportGroupType
      ];
    setState({
      ...state,
      data: data,
    });
  }, [groupDataStore.reportGroupSelectTree, configStore.selectedReference]);

  const treeSelectOnChange = (_currentNode, selectedNodes) => {
    // console.log('onChange::', _currentNode, selectedNodes);
    groupDataStore.applyPendingChanges({
      selectedReportGroups: selectedNodes.map((node) => {
        return node.value;
      }),
    });
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
      <GroupSearchHeader>
        <GroupSearchTitle>
          Select {groupDataStore.getActiveReportGroupTypePrettyName()}s
        </GroupSearchTitle>
      </GroupSearchHeader>
      {dropdownContainer}
    </GroupSearchContainer>
  );
});

export default GroupSearch;
