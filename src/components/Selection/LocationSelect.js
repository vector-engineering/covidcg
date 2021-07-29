import React, { useMemo, useState, useEffect } from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { getNodeFromPath } from '../../utils/location';

import QuestionButton from '../Buttons/QuestionButton';
import LocationItem from './LocationItem';
import {
  ContainerDiv,
  DropdownHeader,
  Title,
  UnselectButton,
  StyledDropdownTreeSelect,
  SelectedLocationsContainer,
} from './LocationSelect.styles';

const LocationSelect = observer(
  ({ selectedLocationNodes, updateSelectedLocationNodes }) => {
    const { locationDataStore } = useStores();

    const [state, setState] = useState({
      data: locationDataStore.selectTree,
    });

    useEffect(() => {
      setState({ ...state, data: locationDataStore.selectTree });
    }, [locationDataStore.selectTree]);

    const onUnselectAll = (e) => {
      e.preventDefault();
      updateSelectedLocationNodes([]);
      locationDataStore.setSelectedNodes([]);
    };

    const treeSelectOnChange = (_currentNode, selectedNodes) => {
      //console.log('onChange::', currentNode, selectedNodes);

      // Since the tree is rendered in a flat state, we need to get all node
      // children from the original data, via. the node paths
      let selectedNodeObjs = selectedNodes.map((node) => {
        // If the path doesn't exist, then we're at the root node
        if (!Object.prototype.hasOwnProperty.call(node, 'path')) {
          return state.data;
        }

        return getNodeFromPath(state.data, node['path']);
      });
      updateSelectedLocationNodes(selectedNodeObjs);
    };
    const treeSelectOnAction = (node, action) => {
      // console.log('onAction::', action, node);

      if (
        Object.prototype.hasOwnProperty.call(action, 'action') &&
        action.action === 'select-all-children'
      ) {
        const selectedNodeObjs = selectedLocationNodes;

        // If this was called on the "All" node, then it won't have a path
        const parentNode = Object.prototype.hasOwnProperty.call(node, 'path')
          ? getNodeFromPath(state.data, node['path'])
          : state.data;

        // Add each child to the selected nodes
        parentNode.children.forEach((child) => {
          selectedNodeObjs.push(child);
        });
        updateSelectedLocationNodes(selectedNodeObjs);
        // Have to update the store's version as well, since
        // we didn't directly click on these new nodes
        locationDataStore.setSelectedNodes(selectedNodeObjs);
      }
    };
    // const treeSelectOnNodeToggleCurrentNode = (currentNode) => {
    //   console.log('onNodeToggle::', currentNode);
    // };
    const onLocationNodeDeselect = (node_path) => {
      let selectedNodeObjs = selectedLocationNodes.filter(
        (node) => node.path !== node_path
      );
      updateSelectedLocationNodes(selectedNodeObjs);
      // Have to update the store's version as well, since
      // we didn't directly click on these new nodes
      locationDataStore.setSelectedNodes(selectedNodeObjs);
    };

    const dropdownContainer = useMemo(
      () => (
        <StyledDropdownTreeSelect
          mode={'hierarchical'}
          data={state.data}
          className="geo-dropdown-tree-select"
          clearSearchOnChange={false}
          keepTreeOnSearch={true}
          keepChildrenOnSearch={true}
          showPartiallySelected={false}
          showDropdown="always"
          inlineSearchInput={true}
          texts={{
            placeholder: 'Search...',
            noMatches: 'No matches found',
          }}
          onChange={treeSelectOnChange}
          onAction={treeSelectOnAction}
          // onNodeToggle={treeSelectOnNodeToggleCurrentNode}
        />
      ),
      [state.data]
    );

    const selectedLocationItems = [];
    selectedLocationNodes.forEach((node) => {
      selectedLocationItems.push(
        <LocationItem
          key={`selected-location-${node.path}`}
          label={node.label}
          path={node.path}
          onDeselect={onLocationNodeDeselect}
        />
      );
    });

    return (
      <ContainerDiv>
        <DropdownHeader>
          <Title>Selected Locations</Title>
          <QuestionButton
            data-tip={`<p>Click on a checkbox to select all sequences from that location.</p><p>You may select multiple levels of locations, i.e., both "USA" and "California".</p><p>Click on the "↳" button to select all of a location's children, without selecting the location itself. i.e., clicking the "↳" next to "All" will select all 6 continents separately.</p><p>Selected countries will appear as boxes underneath the title. Click on the "x" on the right side of each box to de-select that location.</p><p>You can clear all selected locations by clicking "Unselect All"</p>`}
            data-html="true"
            data-place="right"
            data-for="main-tooltip"
          />
          <UnselectButton
            show={selectedLocationNodes.length > 0}
            onClick={onUnselectAll}
          >
            Unselect All
          </UnselectButton>
        </DropdownHeader>
        <SelectedLocationsContainer>
          {selectedLocationItems}
        </SelectedLocationsContainer>
        {dropdownContainer}
      </ContainerDiv>
    );
  }
);

LocationSelect.defaultProps = {
  data: {},
};

LocationSelect.propTypes = {
  data: PropTypes.oneOfType([PropTypes.object, PropTypes.array]).isRequired,
};

LocationSelect.displayName = 'LocationSelect';

export default LocationSelect;
