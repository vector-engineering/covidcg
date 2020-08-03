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

// https://dowjones.github.io/react-dropdown-tree-select/#/story/hoc-readme
// https://dowjones.github.io/react-dropdown-tree-select/#/story/tree-node-paths-hoc
// utility method to assign object paths.
// this path can be used with lodash.get or similar
// to retrieve a deeply nested object
export function assignObjectPaths(obj, stack) {
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

export function getLocationByNameAndLevel(
  selectTree,
  name,
  level,
  check = false
) {
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

    if (check) {
      validNodes.forEach((node) => {
        node.checked = true;
      });
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

export function getNodeFromPath(rootNode, pathStr) {
  // Get children using the object path
  const pathChunks = pathStr.split('.');
  let nodeObj = rootNode;
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
}

export function deselectAll(selectTree) {
  const traverseAndDeselect = (node) => {
    node.checked = false;
    node.children.forEach((child) => {
      traverseAndDeselect(child);
    });
  };
  traverseAndDeselect(selectTree);
}

// Given locObj (a row from location_map), get the name of the most
// specific level of the object, starting from location --> region
// export function getMostSpecificLocationName(locObj) {
//   if (locObj.location !== -1)
// }
