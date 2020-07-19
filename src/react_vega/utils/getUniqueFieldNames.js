export default function getUniqueFieldNames(objects) {
  const fields = new Set();
  objects.forEach((o) => {
    Object.keys(o).forEach((field) => {
      fields.add(field);
    });
  });

  return fields;
}
