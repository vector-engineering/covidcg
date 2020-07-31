import React from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

const AddToSidepanelCheckbox = observer(({ groupKey }) => {
  const { UIStore } = useStores();

  const isSelected = UIStore.sidebarSelectedGroupKeys.includes(groupKey);

  return (
    <input
      type="checkbox"
      id={groupKey}
      name={`add ${groupKey} to sidepanel`}
      checked={isSelected}
      onChange={() => {
        isSelected
          ? UIStore.onRemoveGroupFromSidebar(groupKey)
          : UIStore.onSelectGroupForSidebar(groupKey);
      }}
    ></input>
  );
});

export default AddToSidepanelCheckbox;
