// Recursively look through children for location IDs
function getLocationIdsFromNode(node) {
  if (node === undefined) {
    return [];
  }

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
  return selectTree;
}

export function selectAll(rootNode) {
  const traverseAndSelect = (node) => {
    node.checked = true;
    node.children.forEach((child) => {
      traverseAndSelect(child);
    });
  };
  traverseAndSelect(rootNode);
}
