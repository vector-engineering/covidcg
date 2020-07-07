import React, { useMemo, useState, useEffect } from 'react';
import PropTypes from 'prop-types';
import DropdownTreeSelect from 'react-dropdown-tree-select';
import styled from 'styled-components';
import _ from 'underscore';
import { toJS } from 'mobx';
import { useStores } from '../stores/connect';
import { asyncStates } from '../stores/uiStore';

const ContainerDiv = styled.div`
  margin-top: 5px;
  padding-top: 8px;

  border-top: 1px solid #aaa;
  display: flex;
  flex-direction: column;
  overflow-y: hidden;

  .location-tree-title {
    margin-left: 15px;
  }
`;

const StyledDropdownTreeSelect = styled(DropdownTreeSelect)`
  margin-top: 3px;
  flex-direction: column;
  display: flex;
  overflow-y: hidden;

  ul.tag-list {
    li:first-child {
      span.placeholder:after {
        content: 'None';
        font-size: 0.9rem;
        font-weight: normal;
        font-style: italic;
      }
    }
  }
  .tag {
    background-color: #ffffff;
    border: 1px solid #ccc;
    padding: 3px 6px;
    border-radius: 3px;
    display: inline-block;
    font-weight: normal;
    &:focus-within {
      background-color: #e9e9e9;
      border-color: #a0a0a0;
    }
  }
  span.placeholder {
    background-color: #ffffff;
    font-size: 0em;
  }
  .tag-remove {
    color: #a0a0a0;
    font-size: 1.25em;
    line-height: 100%;
    cursor: pointer;
    background-color: transparent;
    border: none;
    outline: none;
    &:hover,
    &:focus {
      color: #ff5555;
    }
    &.disabled,
    &.readOnly {
      cursor: not-allowed;
    }
  }

  .node > label {
    cursor: pointer;
    margin-left: 2px;
  }

  .tag-list {
    display: flex;
    padding: 0;
    margin: 0;
    flex-wrap: wrap;
  }

  .tag-item {
    display: inline-block;
    margin: 4px;

    .search {
      border: none;
      border-bottom: 1px solid #ccc;
      outline: none;
    }
    &:last-child {
      margin-right: 4px;
    }
  }

  .node {
    list-style: none;
    white-space: nowrap;
    padding: 4px;
    &.leaf.collapsed {
      display: none;
    }
    &.disabled > * {
      color: gray;
      cursor: not-allowed;
    }
    &.match-in-children.hide {
      .node-label {
        opacity: 0.5;
      }
    }
    &.focused {
      background-color: #f4f4f4;
    }
    .toggle {
      white-space: pre;
      margin-right: 4px;
      cursor: pointer;
    }
  }
  .toggle {
    white-space: pre;
    margin-right: 4px;
    cursor:pointer &:after {
      content: ' ';
    }
    &.collapsed:after {
      content: '+';
    }
    &.expanded:after {
      content: '-';
    }
  }

  .searchModeOn {
    .toggle {
      display: none;
    }
  }

  .checkbox-item,
  .radio-item {
    vertical-align: middle;
    margin: 0 4px 0 0;
    &.simple-select {
      display: none;
    }
  }

  .hide:not(.match-in-children) {
    display: none;
  }

  .dropdown {
    width: 100%;
    display: block;
    position: relative;
    flex-direction: column;
    display: flex;
    overflow: hidden;

    a.dropdown-trigger {
      width: calc(100% - 16px);
      border: none;
      padding: 0px 12px;
      line-height: 20px;
      max-height: 200px;
      display: inline-block;
      overflow: auto;
      min-height: 36px;

      &:focus {
        outline: none;
      }
    }

    .dropdown-content {
      position: relative;
      //width: calc(100% - 10px);
      padding: 4px;
      padding-left: 15px;
      padding-right: 15px;
      background-color: #f8f8f8;
      z-index: 1;
      flex-direction: column;
      display: flex;
      overflow-y: hidden;

      input.search {
        font-size: 1em;
        padding: 5px 8px;
        width: calc(100% - 20px);
        border: 1px solid #aaa;
        border-radius: 3px;
        background-color: #ffffff;
        outline: none;
      }

      ul.root {
        margin-top: 5px;
        padding: 0;
        flex-direction: column;
        display: flex;
        overflow-y: scroll;

        i.toggle {
          font-family: monospace;
          font-size: 1.25em;
          font-style: normal;
          font-weight: 800;
          &:hover {
            color: #888888;
          }
          &:focus {
            outline: none;
          }
        }
        .infinite-scroll-component {
          overflow-x: hidden;
        }
      }
    }
  }

  // Custom node styling
  .fa.fa-info {
    margin-left: 0.5em;
    font-weight: normal;
    font-style: normal;
  }
`;

const DropdownContainer = () => {
  const { covidStore, uiStore } = useStores();
  const [state, setState] = useState({
    data: Object.assign(toJS(covidStore.selectTree), { expanded: true }),
    expanded: [],
  });
  const treeSelectOnChange = (currentNode, selectedNodes) => {
    if (
      uiStore.caseDataState === asyncStates.STARTED ||
      uiStore.aggCaseDataState === asyncStates.STARTED
    ) {
      return;
    }
    covidStore.selectLocations(selectedNodes);
  };
  const treeSelectOnAction = (node, action) => {
    console.log('onAction::', action, node);
  };
  const treeSelectOnNodeToggleCurrentNode = (currentNode) => {
    console.log('onNodeToggle::', currentNode);
  };

  useEffect(() => {
    if (!_.isEqual(covidStore.selectTree, state.data)) {
      setState({
        ...state,
        data: Object.assign(toJS(covidStore.selectTree), { expanded: true }),
      });
    }
  }, [covidStore.selectTree]);

  const dropdownContainer = useMemo(
    () => (
      <StyledDropdownTreeSelect
        data={state.data}
        className="geo-dropdown-tree-select"
        clearSearchOnChange={false}
        keepTreeOnSearch={true}
        keepChildrenOnSearch={true}
        showPartiallySelected={true}
        showDropdown="always"
        inlineSearchInput={true}
        texts={{
          placeholder: 'Search...',
          noMatches: 'No matches found',
        }}
        onChange={treeSelectOnChange}
        onAction={treeSelectOnAction}
        onNodeToggle={treeSelectOnNodeToggleCurrentNode}
      />
    ),
    [state.data]
  );

  return (
    <ContainerDiv>
      <span className="location-tree-title">Selected Locations:</span>
      {dropdownContainer}
    </ContainerDiv>
  );
};

DropdownContainer.defaultProps = {
  data: {},
};

DropdownContainer.propTypes = {
  data: PropTypes.oneOfType([PropTypes.object, PropTypes.array]).isRequired,
};

DropdownContainer.displayName = 'DropdownContainer';

export default DropdownContainer;
