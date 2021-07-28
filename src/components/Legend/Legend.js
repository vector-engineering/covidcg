import React, { useState, useEffect } from 'react';
import styled from 'styled-components';
import { toJS } from 'mobx';
import { observer } from 'mobx-react';
import { debounce } from '../../utils/func';

import { useStores } from '../../stores/connect';
import {
  ASYNC_STATES,
  GROUPS,
  COORDINATE_MODES,
  GROUP_SNV,
  DNA_OR_AA,
  SORT_DIRECTIONS,
} from '../../constants/defs.json';
import { COLUMN_NAMES } from './legendUtils';

import SkeletonElement from '../Common/SkeletonElement';
import TableLegend from './TableLegend';

const TableLegendContainer = styled.div`
  width: 100%;
  height: 100vh;
  overflow-y: hidden;
  overflow-x: hidden;
  border-left: 1px #eaeaea solid;
`;

const comparer =
  ({ sortDirection, sortColumn, groupKey, dnaOrAa, coordinateMode }) =>
  (a, b) => {
    // special sorting for snv group
    // If in SNV mode, then sort by position IF:
    // We're in DNA mode OR
    // We're comparing rows that have the same gene/protein
    if (groupKey === GROUP_SNV) {
      let sameGeneOrProtein =
        a.gene_or_protein === b.gene_or_protein &&
        (coordinateMode === COORDINATE_MODES.COORD_GENE ||
          COORDINATE_MODES.COORD_PROTEIN);
      if (
        sortColumn === COLUMN_NAMES.GROUP &&
        (dnaOrAa === DNA_OR_AA.DNA || sameGeneOrProtein)
      ) {
        if (
          sortDirection === SORT_DIRECTIONS.SORT_ASC ||
          sortDirection === SORT_DIRECTIONS.SORT_NONE
        ) {
          return a.pos - b.pos;
        } else {
          return b.pos - a.pos;
        }
      }
    }
    if (
      sortDirection === SORT_DIRECTIONS.SORT_ASC ||
      sortDirection === SORT_DIRECTIONS.SORT_NONE
    ) {
      return a[sortColumn] > b[sortColumn] ? 1 : -1;
    }
    if (sortDirection === SORT_DIRECTIONS.SORT_DESC) {
      return a[sortColumn] < b[sortColumn] ? 1 : -1;
    }
  };

const Legend = observer(() => {
  const { dataStore, UIStore, configStore, groupDataStore, snpDataStore } =
    useStores();

  const [state, setState] = useState({
    legendItems: [],
    sortColumn: COLUMN_NAMES.COUNTS,
    sortDir: SORT_DIRECTIONS.SORT_DESC,
  });

  const updateHoverGroup = debounce((group) => {
    configStore.updateHoverGroup(group);
  }, 10);

  const onClickColumnHeader = ({ columnName }) => {
    if (state.sortColumn === columnName) {
      if (state.sortDir === SORT_DIRECTIONS.SORT_ASC) {
        setState({ ...state, sortDir: SORT_DIRECTIONS.SORT_DESC });
      } else {
        setState({ ...state, sortDir: SORT_DIRECTIONS.SORT_ASC });
      }
    } else {
      setState({
        ...state,
        sortColumn: columnName,
        sortDir: SORT_DIRECTIONS.SORT_DESC,
      });
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
      configStore.selectedGroups.find(
        (group) => group.group === selectedGroup
      ) !== undefined
    ) {
      newGroups = configStore.selectedGroups.filter(
        (group) => !(group.group === selectedGroup)
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
    let legendItems = toJS(dataStore.groupCounts);

    // groupCounts is structured as:
    // [{ group_id, counts }]
    // Where group_id is either a SNV ID (in SNV mode)
    // or a string representing e.g. a lineage

    // Get some additional data:
    // 1) Group Name (get SNV name from SNV ID if in SNV mode)
    // 2) Color
    // 3) Calculate percent based off counts and total sequences
    // 4) SNV gene/protein (for sorting SNVs - AA SNV mode only)
    // 5) SNV position (for sorting SNVs - SNV mode only)
    if (configStore.groupKey === GROUP_SNV) {
      legendItems.forEach((record) => {
        let snv = snpDataStore.intToSnv(
          configStore.dnaOrAa,
          configStore.coordinateMode,
          record.group_id
        );
        record.group = snv.snp_str;
        record.group_name = snv.name;
        record.color = snv.color;
        record.percent =
          record.counts / dataStore.numSequencesAfterAllFiltering;
        record.pos = snv.pos;
        // If we're in DNA mode, then leave this null
        // Otherwise, get the gene or protein, depending on our AA mode
        record.gene_or_protein =
          configStore.dnaOrAa === DNA_OR_AA.DNA
            ? null
            : configStore.coordinateMode === COORDINATE_MODES.COORD_GENE
            ? snv.gene
            : snv.protein;
      });
    } else {
      legendItems.forEach((record) => {
        // For non-SNV groups, the name is the same as the ID
        record.group = record.group_id;
        record.group_name = record.group_id;
        record.color = groupDataStore.getGroupColor(
          configStore.groupKey,
          record.group_id
        );
        record.percent =
          record.counts / dataStore.numSequencesAfterAllFiltering;
      });
    }

    // Set aside the reference, and remove it from the rows list
    // Also set aside the "Other" group, if it exists
    // Sort the list, then add the reference back to the beginning
    // and add the other group back to the end
    const refItem = legendItems.find(
      (item) => item.group === GROUPS.REFERENCE_GROUP
    );
    const otherItem = legendItems.find(
      (item) => item.group === GROUPS.OTHER_GROUP
    );
    legendItems = legendItems.filter(
      (item) =>
        !(
          item.group === GROUPS.REFERENCE_GROUP ||
          item.group === GROUPS.OTHER_GROUP
        )
    );
    legendItems = legendItems.sort(
      comparer({
        sortColumn: state.sortColumn,
        sortDirection: state.sortDir,
        groupKey: configStore.groupKey,
        dnaOrAa: configStore.dnaOrAa,
        coordinateMode: configStore.coordinateMode,
      })
    );
    if (refItem !== undefined) {
      legendItems.unshift(refItem);
    }
    if (otherItem !== undefined) {
      legendItems.push(otherItem);
    }

    return legendItems;
  };

  useEffect(() => {
    let _arr = [...state.legendItems];
    _arr = _arr.sort(
      comparer({
        sortColumn: state.sortColumn,
        sortDirection: state.sortDir,
        groupKey: configStore.groupKey,
        dnaOrAa: configStore.dnaOrAa,
        coordinateMode: configStore.coordinateMode,
      })
    );
    setState({ ...state, legendItems: _arr });
  }, [state.sortColumn, state.sortDir]);

  useEffect(() => {
    if (UIStore.caseDataState !== ASYNC_STATES.SUCCEEDED) {
      return;
    }

    setState({
      ...state,
      legendItems: getLegendKeys().sort(
        comparer({
          sortColumn: state.sortColumn,
          sortDirection: state.sortDir,
          groupKey: configStore.groupKey,
          dnaOrAa: configStore.dnaOrAa,
          coordinateMode: configStore.coordinateMode,
        })
      ),
    });
  }, [UIStore.caseDataState]);

  if (UIStore.caseDataState === ASYNC_STATES.STARTED) {
    return (
      <div
        style={{
          height: '100%',
        }}
      >
        {Array.from({ length: 50 }, (_, i) => (
          <SkeletonElement
            key={`legend-loading-bar-${i}`}
            delay={5 + i + (i % 2) * 12.5}
            height={25}
          />
        ))}
      </div>
    );
  }

  return (
    <TableLegendContainer>
      <TableLegend
        legendItems={state.legendItems}
        updateHoverGroup={updateHoverGroup}
        updateSelectGroup={onItemSelect}
        sortColumn={state.sortColumn}
        sortDir={state.sortDir}
        onClickColumnHeader={onClickColumnHeader}
      />
    </TableLegendContainer>
  );
});

export default Legend;
