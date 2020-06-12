import locations from '../../data/location_map.json';
import selectTree from '../../data/geo_select_tree.json';
// import _ from 'underscore';

// Create ID -> location object/hashmap
const locationIdToNameMap = {};
locations.forEach((loc) => {
  locationIdToNameMap[loc['index']] = loc;
});

export function loadSelectTree() {
  return selectTree;
}

export function getLocationIds(selectedNodes) {
  // Get all locations matching the selected nodes
  let locationIds = [];
  let selectedNode = {};
  let location = {};

  // For each selected node
  for (let i = 0; i < selectedNodes.length; i++) {
    selectedNode = selectedNodes[i];
    // For each location in the cases dataframe
    for (let j = 0; j < locations.length; j++) {
      location = locations[j];

      if (
        selectedNode['level'] === 'region' &&
        selectedNode['value'] === location['region']
      ) {
        locationIds.push(location['index']);
      } else if (
        selectedNode['level'] === 'country' &&
        selectedNode['region'] === location['region'] &&
        selectedNode['value'] === location['country']
      ) {
        locationIds.push(location['index']);
      } else if (
        selectedNode['level'] === 'division' &&
        selectedNode['region'] === location['region'] &&
        selectedNode['country'] === location['country'] &&
        selectedNode['value'] === location['division']
      ) {
        locationIds.push(location['index']);
      } else if (
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
    for (let i = 0; i < nodes.length; i++) {
      let node = nodes[i];
      if (node.level === level && node.value === name) {
        validNodes.push(node);
      }

      // Iterate over all children of child
      validNodes = validNodes.concat(traverseTree(node.children));
    }

    return validNodes;
  }

  return traverseTree(selectTree.children, name, level);
}

export function getLocationNameByIds(locationIds) {
  const out = [];
  locationIds.forEach((locId) => {
    if (typeof locId !== 'string') {
      locId = locId.toString();
    }

    out.push(locationIdToNameMap[locId]);
  });
  return out;
}

// Given locObj (a row from location_map), get the name of the most
// specific level of the object, starting from location --> region
// export function getMostSpecificLocationName(locObj) {
//   if (locObj.location !== -1)
// }
