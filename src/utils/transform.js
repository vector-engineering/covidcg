/*
 * Vega-transform-like functions
 * These are probably really inefficient and could use some
 * performance improvements
 */

// GROUPBY-AGGREGATE
// Tries to be the same as the vega transform of the same name
// groupby = unique keys to group by (one row for each unique combination)
// fields = fields to aggregate
// ops = operations to aggregate by
// as = field to put aggregate value in
export function aggregate({ data, groupby, fields, ops, as }) {
  // Get unique key combinations, and store them in this
  // hierarchical object of objects, one level per group
  // Since this is in a map, populating this should be pretty fast
  const aggContainer = {};

  // Helper to get the current object, depending how deep we are in
  // groupby
  let curContainer;
  let i, j;
  // For each grouping
  for (i = 0; i < groupby.length; i++) {
    data.forEach((row) => {
      // Get the current object corresponding to this row and group level
      // If on the first level, then this is just the master container, agg
      // i.e., if grouped by location and date, and currently on date level,
      //       then this would be agg[location]
      curContainer = aggContainer;
      for (j = 0; j < i; j++) {
        curContainer = curContainer[row[groupby[j]]];
      }

      // Add a key to this object, if it doesn't already exist
      if (!(row[groupby[i]] in curContainer)) {
        curContainer[row[groupby[i]]] = {};
        // If we're not at the final group level, then keep going
        if (i < groupby.length - 1) {
          return;
        }
      }

      // If we're at the final group level, add all the fields we
      // want to aggregate
      if (i === groupby.length - 1) {
        fields.forEach((field) => {
          // Create an array for this field, if the field doesn't exist yet
          if (!(field in curContainer[row[groupby[i]]])) {
            curContainer[row[groupby[i]]][field] = [];
          }
          // Append the current value to this array. We'll do an aggregate
          // function over this list later, depending on whats in the op param
          curContainer[row[groupby[i]]][field].push(row[field]);
        });
      }
    });
  }

  // Use recursion to dig through this hierarchical object,
  // calculate aggregate values, and return a list of objects
  // instead of a nested, hierarchical object
  const recurseThroughKeys = (obj, keys = [], level = 0) => {
    // Defining helper variables first
    let k, _keys, row;
    let rows = [];

    // Get all the keys under the current one
    let groupKeys = Object.keys(obj);

    // For each key
    groupKeys.forEach((key) => {
      // Add the current key to the key chain
      _keys = keys.slice();
      _keys.push(key);

      // If we're on the last level, then calculate aggregate values
      // and return a row with the key chain and aggregate values in
      // the fields as defined in the as param
      if (level === groupby.length - 1) {
        // Create the object for this row to return
        row = {};
        // Add the keys from the key chain
        for (k = 0; k < groupby.length; k++) {
          row[groupby[k]] = _keys[k];
        }
        // Calculate aggregate values for each field
        for (k = 0; k < fields.length; k++) {
          if (ops[k] === 'sum') {
            row[as[k]] = obj[key][fields[k]].reduce(
              (memo, num) => memo + num,
              0
            );
          } else if (ops[k] === 'mean') {
            row[as[k]] =
              obj[key][fields[k]].reduce((memo, num) => memo + num, 0) /
              obj[key][fields[k]].length;
          } else if (ops[k] === 'max') {
            row[as[k]] = obj[key][fields[k]].reduce(
              (memo, num) => (memo = num > memo ? num : memo)
            );
          } else if (ops[k] === 'min') {
            row[as[k]] = obj[key][fields[k]].reduce(
              (memo, num) => (memo = num < memo ? num : memo)
            );
          } else if (ops[k] === 'first') {
            row[as[k]] = obj[key][fields[k]][0];
          }
        }
        // Push this row to the output
        rows.push(row);
      } else {
        // We're not at the last level, keep going.
        // Extend the key chain and increment the level counter
        rows = rows.concat(recurseThroughKeys(obj[key], _keys, level + 1));
      }
    });
    return rows;
  };

  return recurseThroughKeys(aggContainer);
}
