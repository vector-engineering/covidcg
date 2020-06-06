import React from 'react';

const AddToSidepanelCheckbox = ({ id }) => {
  return (
    <input
      type="checkbox"
      id={id}
      name={`add ${id} to sidepanel;`}
      value={`${id} in sidepanel`}
    ></input>
  );
};

export default AddToSidepanelCheckbox;
