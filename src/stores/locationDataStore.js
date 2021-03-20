import { observable, action } from 'mobx';
import { getNodeFromPath, deselectAll } from '../utils/location';
import { asyncDataStoreInstance } from '../components/App';

// https://dowjones.github.io/react-dropdown-tree-select/#/story/hoc-readme
// https://dowjones.github.io/react-dropdown-tree-select/#/story/tree-node-paths-hoc
// utility method to assign object paths.
// this path can be used with lodash.get or similar
// to retrieve a deeply nested object
// function assignObjectPaths(obj, stack) {
//   const isArray = Array.isArray(obj);
//   Object.keys(obj).forEach((k) => {
//     const node = obj[k];
//     const key = isArray ? `[${k}]` : k;

//     if (typeof node === 'object') {
//       node.path = stack ? `${stack}.${key}` : key;
//       assignObjectPaths(node, node.path);
//     }
//   });
// }

// Simplified version, only operate on the "children" field
// instead of on all node keys
function assignObjectPaths(node, stack) {
  if (Object.prototype.hasOwnProperty.call(node, 'children')) {
    node.children.forEach((child, index) => {
      child.path = `${stack ? stack + '.' : ''}children.[${index}]`;
      assignObjectPaths(child, child.path);
    });
  }
}

function addSelectAllChildrenActions(node) {
  if (!Object.prototype.hasOwnProperty.call(node, 'actions')) {
    node.actions = [];
  }

  if (node.children.length > 0) {
    node.actions.push({
      className: 'select-all-children',
      text: '',
      title: 'Select all children',
      action: 'select-all-children',
    });

    node.children.forEach((child) => {
      addSelectAllChildrenActions(child);
    });
  }
}

function recursiveMapIdToStr(map, node) {
  // Add self
  if (Object.prototype.hasOwnProperty.call(node, 'location_id')) {
    if (node.level === 'region') {
      map[node.location_id] = [node.value, '', '', ''];
    } else if (node.level === 'country') {
      map[node.location_id] = [node.region, node.value, '', ''];
    } else if (node.level === 'division') {
      map[node.location_id] = [node.region, node.country, node.value, ''];
    } else if (node.level === 'location') {
      map[node.location_id] = [
        node.region,
        node.country,
        node.division,
        node.value,
      ];
    }
  }
  // Add children
  node.children.forEach((child) => {
    recursiveMapIdToStr(map, child);
  });
}

export class LocationDataStore {
  @observable selectTree = {};
  locationIdToStrMap = {};

  constructor() {}

  init() {
    const selectTree = asyncDataStoreInstance.data.geo_select_tree;
    // Create object paths on each node for easier traversal back and forth
    assignObjectPaths(selectTree);

    // By default, show the tree as expanded so it doesn't only show the "All" node
    selectTree.expanded = true;
    // Also expand the North America node by default
    // selectTree.children.forEach((child) => {
    //   if (child.value === 'North America') {
    //     child.expanded = true;
    //   }
    // });
    addSelectAllChildrenActions(selectTree);

    this.selectTree = selectTree;
    recursiveMapIdToStr(this.locationIdToStrMap, this.selectTree);
  }

  @action
  setSelectedNodes(selectedNodes) {
    // Get the current tree
    const selectTree = Object.assign({}, this.selectTree);
    deselectAll(selectTree);

    // Select the specified nodes
    selectedNodes.forEach((node) => {
      const nodeObj =
        'path' in node ? getNodeFromPath(selectTree, node['path']) : selectTree;
      nodeObj.checked = true;
      // Select all of the nodes children
      // selectAll(nodeObj);
    });

    this.selectTree = selectTree;
  }

  @action
  deselectAll() {
    this.selectTree = deselectAll(Object.assign({}, this.selectTree));
  }

  @action
  getLocationStrFromId(locationId) {
    return this.locationIdToStrMap[locationId];
  }
}
