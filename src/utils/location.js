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

// Recursively look through children for location IDs
function getLocationIdsFromNode(node) {
  return [node.location_id].concat(
    node.children.reduce((memo, child) => {
      return memo.concat(getLocationIdsFromNode(child));
    }, [])
  );
}

export function getLocationIds(selectedNodes) {
  // Get all locations matching the selected nodes
  let locationIds = [];

  // For each selected node
  selectedNodes.forEach((node) => {
    locationIds.push(Array.from(new Set(getLocationIdsFromNode(node))));
  });

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
