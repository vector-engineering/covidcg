import React, { useMemo, useState } from 'react';
import PropTypes from 'prop-types';
import DropdownTreeSelect from 'react-dropdown-tree-select';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';
import { ASYNC_STATES } from '../../constants/UI';

import { loadSelectTree } from '../../utils/location';

const ContainerDiv = styled.div`
  margin-top: 2px;
  padding-top: 8px;

  border-top: 1px solid #aaa;
  display: flex;
  flex-direction: column;
  // overflow-y: hidden;
  overflow-y: scroll;

  .location-tree-title {
    margin-left: 15px;
  }
`;

const StyledDropdownTreeSelect = styled(DropdownTreeSelect)`
  margin-top: 3px;
  flex-direction: column;
  display: flex;
  // overflow-y: hidden;
  overflow-y: scroll;

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
    // overflow: hidden;
    align-items: stretch;

    a.dropdown-trigger {
      width: calc(100% - 16px);
      border: none;
      padding: 0px 12px;
      line-height: 20px;

      &:focus {
        outline: none;
      }
    }

    .dropdown-content {
      position: relative;
      //width: calc(100% - 10px);
      flex-grow: 1;
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

// https://dowjones.github.io/react-dropdown-tree-select/#/story/hoc-readme
// https://dowjones.github.io/react-dropdown-tree-select/#/story/tree-node-paths-hoc
// utility method to assign object paths.
// this path can be used with lodash.get or similar
// to retrieve a deeply nested object
const assignObjectPaths = (obj, stack) => {
  const isArray = Array.isArray(obj);
  Object.keys(obj).forEach((k) => {
    const node = obj[k];
    const key = isArray ? `[${k}]` : k;

    if (typeof node === 'object') {
      node.path = stack ? `${stack}.${key}` : key;
      assignObjectPaths(node, node.path);
    }
  });
};

let initialData = Object.assign(loadSelectTree(), {
  expanded: true,
});
assignObjectPaths(initialData);

const DropdownContainer = () => {
  const { UIStore, configStore } = useStores();

  const [state] = useState({
    data: initialData,
    expanded: [],
  });
  const treeSelectOnChange = (currentNode, selectedNodes) => {
    // Since the tree is rendered in a flat state, we need to get all node
    // children from the original data, via. the node paths
    let selectedNodeObjs = selectedNodes.map((node) => {
      // If the path doesn't exist, then we're at the root node
      if (!Object.prototype.hasOwnProperty.call(node, 'path')) {
        return state.data;
      }

      // Get children using the object path
      const pathChunks = node['path'].split('.');
      let nodeObj = state.data;
      pathChunks.forEach((chunk) => {
        // If this chunk is an array index:
        if (chunk.charAt(0) === '[' && chunk.charAt(chunk.length - 1) === ']') {
          const arrRx = /\[([0-9]*)\]/g;
          const arrIndexMatch = arrRx.exec(chunk);
          nodeObj = nodeObj[parseInt(arrIndexMatch[1])];
        } else {
          nodeObj = nodeObj[chunk];
        }
      });
      return nodeObj;
    });

    if (
      UIStore.caseDataState === ASYNC_STATES.STARTED ||
      UIStore.aggCaseDataState === ASYNC_STATES.STARTED
    ) {
      return;
    }
    configStore.selectLocations(selectedNodeObjs);
  };
  // const treeSelectOnAction = (node, action) => {
  //   console.log('onAction::', action, node);
  // };
  // const treeSelectOnNodeToggleCurrentNode = (currentNode) => {
  //   console.log('onNodeToggle::', currentNode);
  // };

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
        // onAction={treeSelectOnAction}
        // onNodeToggle={treeSelectOnNodeToggleCurrentNode}
      />
    ),
    [state.data]
  );

  return (
    <ContainerDiv>
      <span className="location-tree-title">Selected Locations</span>
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
