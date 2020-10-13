import { observable, toJS, action } from 'mobx';

import { getNodeFromPath, deselectAll, selectAll } from '../utils/location';

import { asyncDataStoreInstance } from '../components/App';

// https://dowjones.github.io/react-dropdown-tree-select/#/story/hoc-readme
// https://dowjones.github.io/react-dropdown-tree-select/#/story/tree-node-paths-hoc
// utility method to assign object paths.
// this path can be used with lodash.get or similar
// to retrieve a deeply nested object
function assignObjectPaths(obj, stack) {
  const isArray = Array.isArray(obj);
  Object.keys(obj).forEach((k) => {
    const node = obj[k];
    const key = isArray ? `[${k}]` : k;

    if (typeof node === 'object') {
      node.path = stack ? `${stack}.${key}` : key;
      assignObjectPaths(node, node.path);
    }
  });
}

export class LocationDataStore {
  locations = [];
  @observable selectTree = {};

  // ID -> location object/hashmap
  locationIdToNameMap = {};

  constructor() {}

  init() {
    this.locations = asyncDataStoreInstance.data.location_map;
    const selectTree = asyncDataStoreInstance.data.geo_select_tree;
    // Create object paths on each node for easier traversal back and forth
    assignObjectPaths(selectTree);
    // By default, show the tree as expanded so it doesn't only show the "All" node
    selectTree.expanded = true;
    this.selectTree = selectTree;

    // Create ID -> location object/hashmap
    this.locations.forEach((loc) => {
      this.locationIdToNameMap[loc['index']] = loc;
    });
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
      // Select all of the nodes children
      selectAll(nodeObj);
    });

    this.selectTree = selectTree;
  }

  @action
  deselectAll() {
    this.selectTree = deselectAll(toJS(this.selectTree));
  }
}
