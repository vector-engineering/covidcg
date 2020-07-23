import React from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

const AddToSidepanelCheckbox = observer(({ groupKey }) => {
  const { uiStore } = useStores();

  const isSelected = uiStore.sidebarSelectedGroupKeys.includes(groupKey);

  return (
    <input
      type="checkbox"
      id={groupKey}
      name={`add ${groupKey} to sidepanel`}
      checked={isSelected}
      onChange={() => {
        isSelected
          ? uiStore.onRemoveGroupFromSidebar(groupKey)
          : uiStore.onSelectGroupForSidebar(groupKey);
      }}
    ></input>
  );
});

export default AddToSidepanelCheckbox;
