import React, { useMemo, useState, useEffect } from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { getNodeFromPath } from '../../utils/location';

import DeselectButton from '../Buttons/DeselectButton';
import QuestionButton from '../Buttons/QuestionButton';

import {
  ContainerDiv,
  DropdownHeader,
  Title,
  UnselectButton,
  SelectedLocationsContainer,
  LocationItemContainer,
  LocationItemLabel,
  StyledDropdownTreeSelect,
} from './LocationSelect.styles';

const LocationItem = ({ label, path, onDeselect }) => {
  const handleDeselect = (e) => {
    e.preventDefault();
    onDeselect(path);
  };

  return (
    <LocationItemContainer>
      <LocationItemLabel>{label}</LocationItemLabel>
      <DeselectButton title={'Deselect'} onClick={handleDeselect}>
        ×
      </DeselectButton>
    </LocationItemContainer>
  );
};
LocationItem.propTypes = {
  label: PropTypes.string.isRequired,
  path: PropTypes.string,
  onDeselect: PropTypes.func.isRequired,
};

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
      let allRegions = false;

      // Since the tree is rendered in a flat state, we need to get all node
      // children from the original data, via. the node paths
      let selectedNodeObjs = selectedNodes.map((node) => {
        // If the path doesn't exist, then we're at the root node
        if (!Object.prototype.hasOwnProperty.call(node, 'path')) {
          // Select all regions when All is checked
          allRegions = true;
          return state.data;
        }

        return getNodeFromPath(state.data, node['path']);
      });

      if (allRegions) {
        // Add any regions not already in selectedNodeObjs, and remove self
        state.data.children.forEach((region) => {
          if (!selectedNodeObjs.includes(region)) {
            selectedNodeObjs.push(region);
          }
        });
        // Remove self
        selectedNodeObjs = selectedNodeObjs.filter((node) => {
          return node.value != 'All';
        });
      }

      updateSelectedLocationNodes(selectedNodeObjs);
      // Have to update the store's version as well, since
      // we didn't directly click on these new nodes
      locationDataStore.setSelectedNodes(selectedNodeObjs);
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

    // Maintain tree expansion state
    const treeSelectOnNodeToggleCurrentNode = (currentNode) => {
      // console.log(currentNode);
      const data = Object.assign({}, state.data);

      // Recursively go through and find the node to expand
      const traverseAndExpand = (node) => {
        if (node.path === currentNode.path) {
          node.expanded = currentNode.expanded;
        }

        if ('children' in node) {
          node.children.forEach((child) => {
            traverseAndExpand(child);
          });
        }
      };
      traverseAndExpand(data);

      setState({
        ...state,
        data,
      });
    };

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
          onNodeToggle={treeSelectOnNodeToggleCurrentNode}
        />
      ),
      [state.data]
    );

    const selectedLocationItems = [];
    selectedLocationNodes.forEach((node) => {
      if (!selectedLocationItems.includes(node)) {
        selectedLocationItems.push(
          <LocationItem
            key={`selected-location-${node.path}`}
            label={node.label}
            path={node.path}
            onDeselect={onLocationNodeDeselect}
          />
        );
      }
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
