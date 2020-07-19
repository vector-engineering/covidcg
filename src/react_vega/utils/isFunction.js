export default function isFunction(functionToCheck) {
  const getType = {};
  return (
    !!functionToCheck &&
    getType.toString.call(functionToCheck) === '[object Function]'
  );
}
