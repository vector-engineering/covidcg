import React, { useState, useEffect } from 'react';
import styled from 'styled-components';
import { observer } from 'mobx-react';
import _ from 'underscore';

import { useStores } from '../../stores/connect';
import {
  ASYNC_STATES,
  GROUPS,
  COORDINATE_MODES,
  GROUP_SNV,
  DNA_OR_AA,
} from '../../constants/defs.json';
import { sortLegendItems } from './legendutils';
import TableLegend from './TableLegend';

const TableLegendContainer = styled.div`
  width: 100%;
  height: 100%;
`;

const comparer = ({
  sortDirection,
  sortColumn,
  groupKey,
  dnaOrAa,
  coordinateMode,
}) => (a, b) => {
  // special sorting for snv group
  // If in SNV mode, then sort by position IF:
  // We're in DNA mode OR
  // We're comparing rows that have the same gene/protein
  if (groupKey === GROUP_SNV) {
    let sameGeneOrProtein =
      (a.gene === b.gene && coordinateMode === COORDINATE_MODES.COORD_GENE) ||
      (a.protein === b.protein &&
        coordinateMode === COORDINATE_MODES.COORD_PROTEIN);

    if (
      sortColumn === 'group' &&
      (dnaOrAa === DNA_OR_AA.DNA || sameGeneOrProtein) &&
      sortDirection === 'ASC'
    ) {
      return a.pos - b.pos;
    }
    if (
      sortColumn === 'group' &&
      (dnaOrAa === DNA_OR_AA.DNA || sameGeneOrProtein) &&
      sortDirection === 'DESC'
    ) {
      return b.pos - a.pos;
    }
  }
  if (sortDirection === 'ASC' || sortDirection === 'None') {
    return a[sortColumn] > b[sortColumn] ? 1 : -1;
  }
  if (sortDirection === 'DESC') {
    return a[sortColumn] < b[sortColumn] ? 1 : -1;
  }
};

const LegendContainer = observer(() => {
  const { dataStore, UIStore, configStore } = useStores();

  const [legendItems, setLegendItems] = useState([]);
  const [sortColumn, setSortColumn] = useState('cases_percent');
  const [sortDir, setSortDir] = useState('DESC');

  const updateHoverGroup = _.debounce((group) => {
    configStore.updateHoverGroup(group);
  }, 10);

  const onClickColumnHeader = ({ columnName }) => {
    if (sortColumn === columnName && sortDir === 'DESC') {
      setSortDir('ASC');
    } else {
      setSortColumn(columnName);
      setSortDir('DESC');
    }
  };

  const onItemSelect = (e) => {
    const selectedGroup = e.target.getAttribute('data-group');
    console.log(selectedGroup);
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
    const refItem = _.findWhere(_legendItems, {
      group: GROUPS.REFERENCE_GROUP,
    });
    const otherItem = _.findWhere(_legendItems, { group: GROUPS.OTHER_GROUP });
    _legendItems = _.reject(
      _legendItems,
      (item) =>
        item.group === GROUPS.REFERENCE_GROUP ||
        item.group === GROUPS.OTHER_GROUP
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
    let _arr = [...legendItems];
    _arr = _arr.sort(
      comparer({
        sortColumn,
        sortDirection: sortDir,
        groupKey: configStore.groupKey,
        dnaOrAa: configStore.dnaOrAa,
        coordinateMode: configStore.coordinateMode,
      })
    );
    setLegendItems(_arr);
  }, [sortColumn, sortDir]);

  useEffect(() => {
    if (UIStore.caseDataState !== ASYNC_STATES.SUCCEEDED) {
      return;
    }

    setLegendItems(
      getLegendKeys().sort(
        comparer({
          sortColumn,
          sortDirection: sortDir,
          groupKey: configStore.groupKey,
          dnaOrAa: configStore.dnaOrAa,
          coordinateMode: configStore.coordinateMode,
        })
      )
    );
  }, [UIStore.caseDataState]);

  if (UIStore.caseDataState !== ASYNC_STATES.SUCCEEDED) {
    return <div />;
  }

  return (
    <TableLegendContainer>
      <TableLegend
        legendItems={legendItems}
        updateHoverGroup={updateHoverGroup}
        updateSelectGroup={onItemSelect}
        sortColumn={sortColumn}
        sortDir={sortDir}
        onClickColumnHeader={onClickColumnHeader}
      />
    </TableLegendContainer>
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
