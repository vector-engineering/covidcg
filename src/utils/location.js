import locations from '../../processed_data/locations.json';
import select_tree from '../../processed_data/geo_select_tree.json';

export function loadSelectTree() {
  return select_tree;
}

export function getLocationIds(selectedNodes) {
  // Get all locations matching the selected nodes
  let locationIds = [];
  let selectedNode = {};
  let location = {};

  // For each selected node
  for(let i = 0; i < selectedNodes.length; i++) {
    selectedNode = selectedNodes[i];
    // For each location in the cases dataframe
    for(let j = 0; j < locations.length; j++) {
      location = locations[j];

      if(
        selectedNode['level'] === 'region' && 
        selectedNode['value'] === location['region']
      ) {
        locationIds.push(location['index']);
      }
      else if(
        selectedNode['level'] === 'country' && 
        selectedNode['region'] === location['region'] &&
        selectedNode['value'] === location['country']
      ) {
        locationIds.push(location['index']);
      }
      else if(
        selectedNode['level'] === 'division' && 
        selectedNode['region'] === location['region'] &&
        selectedNode['country'] === location['country'] &&
        selectedNode['value'] === location['division']
      ) {
        locationIds.push(location['index']);
      }
      else if(
        selectedNode['level'] === 'location' && 
        selectedNode['region'] === location['region'] &&
        selectedNode['country'] === location['country'] &&
        selectedNode['division'] === location['division'] &&
        selectedNode['value'] === location['location']
      ) {
        locationIds.push(location['index']);
      }
    }
  }

  return locationIds;
}

export function getLocationByNameAndLevel(selectTree, name, level) {
  
  // Traverse tree with recursion
  function traverseTree(nodes) {
    let validNodes = [];

    // Check all children
    for(let i = 0; i < nodes.length; i++) {
      let node = nodes[i];
      if(node.level === level && node.value === name) {
        validNodes.push(node)
      }

      // Iterate over all children of child
      validNodes = validNodes.concat(traverseTree(node.children));
    }

    return validNodes;
  }

  return traverseTree(selectTree.children, name, level);
}
