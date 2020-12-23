import React, { useState, useEffect } from 'react';
import { observer } from 'mobx-react';
import _ from 'underscore';

import { useStores } from '../../stores/connect';
import { ASYNC_STATES } from '../../constants/UI';
import { REFERENCE_GROUP, OTHER_GROUP } from '../../constants/groups';
import { sortLegendItems } from './legendutils';
import TableLegend from './TableLegend';

const LegendContainer = observer(() => {
  const { dataStore, UIStore, configStore } = useStores();

  const [legendItems, setLegendItems] = useState([]);

  const updateHoverGroup = _.debounce((group) => {
    configStore.updateHoverGroup(group);
  }, 10);

  const onItemSelect = (e) => {
    const selectedGroup = e.target.getAttribute('data-group');
    let newGroups;

    // If the click was not on an item, then unset the selection
    if (selectedGroup === null) {
      newGroups = [];
    }
    // If the item is already selected, then deselect it
    else if (
      _.findWhere(configStore.selectedGroups, { group: selectedGroup }) !==
      undefined
    ) {
      newGroups = _.reject(
        configStore.selectedGroups,
        (group) => group.group == selectedGroup
      );
    } else {
      // Otherwise, add it
      newGroups = [{ group: selectedGroup }];
      // If shift is pressed, then add it to the existing selected groups
      if (UIStore.isKeyPressed(16)) {
        newGroups = newGroups.concat(configStore.selectedGroups);
      }
    }

    configStore.updateSelectedGroups(newGroups);
  };

  const getLegendKeys = () => {
    // Make own copy of the elements, and sort by group
    let _legendItems = JSON.parse(JSON.stringify(dataStore.dataAggGroup));

    // Set aside the reference, and remove it from the rows list
    // Also set aside the "Other" group, if it exists
    // Sort the list, then add the reference back to the beginning
    // and add the other group back to the end
    const refItem = _.findWhere(_legendItems, { group: REFERENCE_GROUP });
    const otherItem = _.findWhere(_legendItems, { group: OTHER_GROUP });
    _legendItems = _.reject(
      _legendItems,
      (item) => item.group === REFERENCE_GROUP || item.group === OTHER_GROUP
    );
    _legendItems = _legendItems.sort(
      sortLegendItems.bind(
        this,
        configStore.groupKey,
        configStore.dnaOrAa,
        configStore.coordinateMode
      )
    );
    if (refItem !== undefined) {
      _legendItems.unshift(refItem);
    }
    if (otherItem !== undefined) {
      _legendItems.push(otherItem);
    }

    return _legendItems;
  };

  useEffect(() => {
    if (UIStore.caseDataState !== ASYNC_STATES.SUCCEEDED) {
      return;
    }

    setLegendItems(getLegendKeys());
  }, [UIStore.caseDataState]);

  console.log(legendItems);

  return (
    <TableLegend
      legendItems={legendItems}
      updateHoverGroup={updateHoverGroup}
      updateSelectGroup={onItemSelect}
    />
  );

  /*
  return <VegaLegend
      legendItems={legendItems}
      updateHoverGroup={updateHoverGroup}
      updateSelectGroup={onItemSelect} 
      />;
  */
});

export default LegendContainer;
