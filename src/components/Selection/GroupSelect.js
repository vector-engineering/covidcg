import React, { useMemo, useState, useEffect } from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { config } from '../../config';

import QuestionButton from '../Buttons/QuestionButton';
import DeselectButton from '../Buttons/DeselectButton';
import StyledDropdownTreeSelect from '../Common/StyledDropdownTreeSelect';

import {
  GroupSelectContainer,
  GroupSelectForm,
  GroupSelectedItemsContainer,
  GroupTitle,
  GroupSelectedItemList,
  GroupSelectedItemContainer,
  GroupSelectedItemLabel,
} from './GroupSelect.styles';

const sortGroupNames = (a, b) => {
  if (a.name < b.name) {
    return -1;
  }
  if (a.name > b.name) {
    return 1;
  }
  return 0;
};

const GroupSelectedItem = ({ group, onDeselect }) => {
  const handleDeselect = (e) => {
    e.preventDefault();
    onDeselect(group);
  };

  return (
    <GroupSelectedItemContainer>
      <GroupSelectedItemLabel>{group}</GroupSelectedItemLabel>
      <DeselectButton title={'Deselect'} onClick={handleDeselect}>
        ×
      </DeselectButton>
    </GroupSelectedItemContainer>
  );
};
GroupSelectedItem.propTypes = {
  group: PropTypes.string.isRequired,
  onDeselect: PropTypes.func.isRequired,
};

const GroupSelect = observer(
  ({ selectedGroupFields, updateSelectedGroupFields }) => {
    // console.log(selectedGroupFields);

    const { groupDataStore } = useStores();

    const createGroupSelectTree = (groups) => {
      const reportGroupSelectTree = {};

      // Construct selection trees
      Object.keys(groups).forEach((groupName) => {
        reportGroupSelectTree[groupName] = [];
        groups[groupName].sort(sortGroupNames).forEach((group) => {
          reportGroupSelectTree[groupName].push({
            label: group.name,
            value: group.name,
            checked:
              Object.prototype.hasOwnProperty.call(
                selectedGroupFields,
                groupName
              ) && selectedGroupFields[groupName].includes(group.name),
          });
        });
      });

      return reportGroupSelectTree;
    };

    const [state, setState] = useState({
      data: createGroupSelectTree(groupDataStore.groups),
    });

    useEffect(() => {
      setState({
        ...state,
        data: createGroupSelectTree(groupDataStore.groups),
      });
    }, [groupDataStore.groups]);

    const treeSelectOnChange = (groupKey, _currentNode, selectedNodes) => {
      // console.log('onChange::', groupKey, _currentNode, selectedNodes);

      const newSelectedGroupFields = JSON.parse(
        JSON.stringify(selectedGroupFields)
      );
      newSelectedGroupFields[groupKey] = selectedNodes.map(
        (node) => node.value
      );

      updateSelectedGroupFields(newSelectedGroupFields);

      // Manually update the tree
      const treeData = JSON.parse(JSON.stringify(state.data));
      treeData[groupKey].forEach((node) => {
        if (newSelectedGroupFields[groupKey].includes(node.value)) {
          node.checked = true;
        }
      });
      setState({ ...state, data: treeData });
    };

    const onNodeDeselect = (groupKey, group) => {
      // console.log(groupKey, group);
      const newSelectedGroupFields = JSON.parse(
        JSON.stringify(selectedGroupFields)
      );
      newSelectedGroupFields[groupKey] = newSelectedGroupFields[
        groupKey
      ].filter((selectedGroup) => selectedGroup !== group);

      updateSelectedGroupFields(newSelectedGroupFields);

      // Manually update the tree
      const treeData = JSON.parse(JSON.stringify(state.data));
      treeData[groupKey].forEach((node) => {
        if (node.value === group) {
          node.checked = false;
        }
      });
      setState({ ...state, data: treeData });
    };

    // const treeSelectOnAction = (node, action) => {
    //   console.log('onAction::', action, node);
    // };

    const groupDropdowns = [];
    Object.keys(config['group_cols']).forEach((groupKey) => {
      groupDropdowns.push(
        useMemo(
          () => (
            <StyledDropdownTreeSelect
              mode={'multiSelect'}
              data={state.data[groupKey]}
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
              onChange={treeSelectOnChange.bind(this, groupKey)}
              // onAction={treeSelectOnAction}
              // onNodeToggle={treeSelectOnNodeToggleCurrentNode}
            />
          ),
          [state.data[groupKey]]
        )
      );
    });

    const groupDialogs = [];
    Object.keys(config['group_cols']).forEach((groupKey, i) => {
      const groupSelectedItems = [];
      if (Object.prototype.hasOwnProperty.call(selectedGroupFields, groupKey)) {
        selectedGroupFields[groupKey].forEach((group) => {
          groupSelectedItems.push(
            <GroupSelectedItem
              key={`selected-group-item-${group}`}
              group={group}
              onDeselect={onNodeDeselect.bind(this, groupKey)}
            />
          );
        });
      }

      groupDialogs.push(
        <GroupSelectForm key={`group-select-item-${groupKey}`}>
          <GroupSelectedItemsContainer>
            <GroupTitle>{config['group_cols'][groupKey].title}:</GroupTitle>
            <GroupSelectedItemList>
              {groupSelectedItems.length > 0 && groupSelectedItems}
              {groupSelectedItems.length === 0 && (
                <span>No {config['group_cols'][groupKey].title}s selected</span>
              )}
            </GroupSelectedItemList>
          </GroupSelectedItemsContainer>
          {groupDropdowns[i]}
        </GroupSelectForm>
      );
    });

    return (
      <GroupSelectContainer>
        <span className="title">
          Selected Groups
          <QuestionButton
            data-tip="<p>Filter for sequences based on their assignment to specified phylogenetic groups</p>"
            data-html={true}
            data-place="left"
            data-for="main-tooltip"
          />
        </span>

        {groupDialogs}
      </GroupSelectContainer>
    );
  }
);

GroupSelect.propTypes = {
  selectedGroupFields: PropTypes.object,
  updateSelectedGroupFields: PropTypes.func,
};

export default GroupSelect;
