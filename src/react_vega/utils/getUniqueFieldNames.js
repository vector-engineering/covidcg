export default function getUniqueFieldNames(objects) {
  const fields = new Set();
  objects.forEach((o) => {
    if (o === undefined || o === null || typeof o !== 'object') {
      return;
    }

    Object.keys(o).forEach((field) => {
      fields.add(field);
    });
  });

  return fields;
}
