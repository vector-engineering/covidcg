import React, { useMemo, useState, useEffect } from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { getNodeFromPath } from '../../utils/location';

import {
  ContainerDiv,
  DropdownHeader,
  UnselectButton,
  StyledDropdownTreeSelect,
} from './LocationSelect.styles';

import { ASYNC_STATES } from '../../constants/defs.json';

const LocationSelect = observer(() => {
  const { UIStore, configStore, locationDataStore } = useStores();

  const [state, setState] = useState({
    data: locationDataStore.selectTree,
  });

  useEffect(() => {
    setState({ ...state, data: locationDataStore.selectTree });
  }, [locationDataStore.selectTree]);

  const onUnselectAll = (e) => {
    e.preventDefault();
    configStore.selectLocations([]);
  };

  const treeSelectOnChange = (currentNode, selectedNodes) => {
    // Since the tree is rendered in a flat state, we need to get all node
    // children from the original data, via. the node paths
    let selectedNodeObjs = selectedNodes.map((node) => {
      // If the path doesn't exist, then we're at the root node
      if (!Object.prototype.hasOwnProperty.call(node, 'path')) {
        return state.data;
      }

      return getNodeFromPath(state.data, node['path']);
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
      <DropdownHeader>
        <span className="location-tree-title">Selected Locations</span>
        <UnselectButton
          show={configStore.selectedLocationNodes.length > 0}
          onClick={onUnselectAll}
        >
          Unselect All
        </UnselectButton>
      </DropdownHeader>
      {dropdownContainer}
    </ContainerDiv>
  );
});

LocationSelect.defaultProps = {
  data: {},
};

LocationSelect.propTypes = {
  data: PropTypes.oneOfType([PropTypes.object, PropTypes.array]).isRequired,
};

LocationSelect.displayName = 'LocationSelect';

export default LocationSelect;
