import React, { useState, useEffect } from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import _ from 'underscore';

import { ASYNC_STATES } from '../../constants/UI';
import {
  GROUP_KEYS,
  DNA_OR_AA,
  COORDINATE_MODES,
} from '../../constants/config';
import { REFERENCE_GROUP, OTHER_GROUP } from '../../constants/groups';

import SkeletonElement from '../Common/SkeletonElement';
import { LegendList } from './Legend.styles';
import LegendItemWrapper from './LegendItem';

const LegendList = observer(() => {
  const { dataStore, UIStore, configStore } = useStores();

  const [state, setState] = useState({
    legendItems: [],
  });
  console.log(state);

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

  const updateHoverGroup = _.debounce((group) => {
    configStore.updateHoverGroup(group);
  }, 10);

  const renderLegendKeys = () => {
    // Make own copy of the elements, and sort by group
    let legendItems = JSON.parse(JSON.stringify(dataStore.dataAggGroup));

    // Set aside the reference, and remove it from the rows list
    // Also set aside the "Other" group, if it exists
    // Sort the list, then add the reference back to the beginning
    // and add the other group back to the end
    const refItem = _.findWhere(legendItems, { group: REFERENCE_GROUP });
    const otherItem = _.findWhere(legendItems, { group: OTHER_GROUP });
    legendItems = _.reject(
      legendItems,
      (item) => item.group === REFERENCE_GROUP || item.group === OTHER_GROUP
    );
    legendItems = legendItems.sort(
      sortLegendItems.bind(
        this,
        configStore.groupKey,
        configStore.dnaOrAa,
        configStore.coordinateMode
      )
    );
    if (refItem !== undefined) {
      legendItems.unshift(refItem);
    }
    if (otherItem !== undefined) {
      legendItems.push(otherItem);
    }

    return legendItems.map((obj) => {
      if (!obj.color) {
        return null;
      }

      return (
        <LegendItemWrapper
          key={`legend-item-${obj.group}`}
          color={obj.color}
          group={obj.group}
          updateHoverGroup={updateHoverGroup}
        />
      );
    });
  };

  const sortLegendItems = (groupKey, dnaOrAa, coordinateMode, a, b) => {
    // If we're grouping by lineage or clade, then sort alphabetically
    // on the lineage/clade
    if (
      groupKey === GROUP_KEYS.GROUP_LINEAGE ||
      groupKey === GROUP_KEYS.GROUP_CLADE
    ) {
      return a.group > b.group ? 1 : -1;
    } else if (groupKey === GROUP_KEYS.GROUP_SNV) {
      // If we're grouping by SNV, figure out whether we're in DNA or AA mode
      if (dnaOrAa === DNA_OR_AA.DNA) {
        // If we're grouping by DNA SNV, then sort by position
        return a.pos > b.pos ? 1 : -1;
      } else {
        // If we're grouping by AA SNV, determine if we're in gene/protein mode
        if (coordinateMode === COORDINATE_MODES.COORD_GENE) {
          // If the same gene, then sort by position
          if (a.gene === b.gene) {
            return a.pos > b.pos ? 1 : -1;
          }
          // Otherwise, sort by gene
          else {
            return a.gene > b.gene ? 1 : -1;
          }
        } else if (coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
          // If the same protein, then sort by position
          if (a.protein === b.protein) {
            return a.pos > b.pos ? 1 : -1;
          }
          // Otherwise, sort by protein
          else {
            return a.protein.toLowerCase() > b.protein.toLowerCase() ? 1 : -1;
          }
        }
      }
    }
  };

  useEffect(() => {
    if (UIStore.caseDataState !== ASYNC_STATES.SUCCEEDED) {
      return;
    }

    setState({ ...state, legendItems: renderLegendKeys() });
  }, [UIStore.caseDataState]);

  // console.log('RE-RENDERING LEGEND');

  if (UIStore.caseDataState === ASYNC_STATES.STARTED) {
    return (
      <div
        style={{
          paddingTop: '12px',
          paddingRight: '24px',
          paddingLeft: '12px',
          paddingBottom: '0px',
        }}
      >
        <SkeletonElement delay={1} height={50} />
      </div>
    );
  }

  return (
    <LegendList onMouseDown={onItemSelect}>{state.legendItems}</LegendList>
  );
});

LegendList.displayName = 'LegendList';

export default LegendList;
