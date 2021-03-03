/* eslint-disable react/display-name */
import React, { useState, useEffect } from 'react';
import styled from 'styled-components';

import EmptyDataTable from './EmptyDataTable';
import 'react-data-grid/dist/react-data-grid.css';
import { observer } from 'mobx-react';
import _ from 'underscore';

import { useStores } from '../../stores/connect';
import { meetsContrastGuidelines } from 'polished';
import { nanmax, nanmin } from '../../utils/math';
import {
  snapGeneNTColors,
  shingAAColors,
  clustalXAAColors,
  zappoAAColors,
  transmembraneAAColors,
} from '../../constants/colors';
import {
  DNA_OR_AA,
  COLOR_MODES,
  COMPARE_COLORS,
  SORT_DIRECTIONS,
  ASYNC_STATES,
  GROUPS,
} from '../../constants/defs.json';
import { config } from '../../config';

import TableOptions from './TableOptions';
import {
  groupColumn,
  getDefaultColumns,
  getSinglePosColumn,
} from './columnDefs';

import SkeletonElement from '../Common/SkeletonElement';
import DataTable from './DataTable';
import RowRenderer from './RowRenderer';

const DataTableContainer = styled.div`
  display: flex;
  flex-direction: column;
  flex-grow: 1;
  position: relative;
  span.position-title {
    font-size: 12px;
    font-weight: 500;
  }
`;

const comparer = ({ sortDirection, sortColumn }) => (a, b) => {
  if (sortDirection === 'ASC' || sortDirection === 'None') {
    return a[sortColumn] > b[sortColumn] ? 1 : -1;
  }
  if (sortDirection === 'DESC') {
    return a[sortColumn] < b[sortColumn] ? 1 : -1;
  }
};

const sortRows = (rows, sortFn) => {
  // Set aside the reference, and remove it from the rows list
  let refRow = _.findWhere(rows, { group: GROUPS.REFERENCE_GROUP });
  rows = _.reject(rows, (row) => row.group == GROUPS.REFERENCE_GROUP);
  rows = rows.sort(sortFn);
  if (refRow !== undefined) {
    rows.unshift(refRow);
  }
  return rows;
};

const NewLineageDataTable = observer(() => {
  const { dataStore, UIStore, configStore, plotSettingsStore } = useStores();

  const buildColumns = () => {
    let _columns = [];
    // For lineage grouping, add lineage column
    _columns.push(
      groupColumn({ title: config.group_cols[configStore.groupKey].title })
    );

    // Get the maximum and minimum counts and percent for the colormaps
    // Ignore those values for the reference row (which are NaN)
    let maxCasesSum = _.reduce(
      dataStore.dataAggGroup,
      (memo, group) => nanmax(memo, group.counts),
      0
    );
    let minCasesSum = _.reduce(
      dataStore.dataAggGroup,
      (memo, group) => nanmin(memo, group.counts),
      0
    );
    let maxCasesPercent = _.reduce(
      dataStore.dataAggGroup,
      (memo, group) => nanmax(memo, group.percent),
      0
    );
    let minCasesPercent = _.reduce(
      dataStore.dataAggGroup,
      (memo, group) => nanmin(memo, group.percent),
      0
    );

    // example row
    // percent: 0.6880290205562273
    // counts: 569
    // group: "B.1"
    // pos_23402: "G"
    // pos_23730: "C"
    // pos_24033: "C"

    _columns = _columns.concat(
      getDefaultColumns({
        minCasesPercent,
        maxCasesPercent,
        minCasesSum,
        maxCasesSum,
      })
    );

    // Build a column for each changing position
    let refRow = _.findWhere(dataStore.dataAggGroup, {
      group: GROUPS.REFERENCE_GROUP,
    });
    if (refRow === undefined) {
      return null;
    }

    Object.keys(refRow).forEach((col) => {
      // Only process columns starting with "pos_"
      if (!col.startsWith('pos_')) {
        return;
      }

      let colors = null;
      let textColors = null;
      if (
        configStore.dnaOrAa === DNA_OR_AA.DNA &&
        (plotSettingsStore.tableCompareColor ===
          COMPARE_COLORS.COLOR_MODE_CODE ||
          plotSettingsStore.tableColorMode === COLOR_MODES.COLOR_MODE_CODE)
      ) {
        colors = snapGeneNTColors;
      } else {
        if (
          plotSettingsStore.tableCompareColor ===
            COMPARE_COLORS.COLOR_MODE_CODE ||
          plotSettingsStore.tableColorMode === COLOR_MODES.COLOR_MODE_CODE
        ) {
          colors = shingAAColors;
        } else if (
          plotSettingsStore.tableCompareColor ===
            COMPARE_COLORS.COLOR_MODE_CLUSTAL ||
          plotSettingsStore.tableColorMode === COLOR_MODES.COLOR_MODE_CLUSTAL
        ) {
          colors = clustalXAAColors;
        } else if (
          plotSettingsStore.tableCompareColor ===
            COMPARE_COLORS.COLOR_MODE_ZAPPO ||
          plotSettingsStore.tableColorMode === COLOR_MODES.COLOR_MODE_ZAPPO
        ) {
          colors = zappoAAColors;
        } else if (
          plotSettingsStore.tableCompareColor ===
            COMPARE_COLORS.COLOR_MODE_ZHAO_LONDON ||
          plotSettingsStore.tableColorMode ===
            COLOR_MODES.COLOR_MODE_ZHAO_LONDON
        ) {
          colors = transmembraneAAColors;
        }
      }

      // Calculate text colorings
      if (colors !== null) {
        textColors = {};
        let scores;
        Object.keys(colors).forEach((key) => {
          scores = meetsContrastGuidelines(colors[key], '#fff');
          textColors[key] = scores['AALarge'] ? '#fff' : '#000';
        });
      }

      // 0-indexed to 1-indexed
      let pos = parseInt(col.substring(4));
      if (configStore.dnaOrAa === DNA_OR_AA.DNA) {
        pos += 1;
      }

      _columns.push(
        getSinglePosColumn({
          pos,
          col,
          colorMode: plotSettingsStore.tableColorMode,
          refRow,
          compareMode: plotSettingsStore.tableCompareMode,
          compareColor: plotSettingsStore.tableCompareColor,
          colors,
          textColors,
        })
      );
    });

    // console.log(_columns);

    return _columns;
  };

  const [state, setState] = useState({
    columns: buildColumns() || [],
    rows: dataStore.dataAggGroup,
  });

  useEffect(() => {
    if (UIStore.caseDataState !== ASYNC_STATES.SUCCEEDED) {
      return;
    }

    setState({
      ...state,
      columns: buildColumns() || [],
      rows: sortRows(
        dataStore.dataAggGroup,
        comparer({
          sortDirection: plotSettingsStore.tableSortDirection,
          sortColumn: plotSettingsStore.tableSortColumn,
        })
      ),
    });
  }, [
    UIStore.caseDataState,
    plotSettingsStore.tableColorMode,
    plotSettingsStore.tableCompareMode,
    plotSettingsStore.tableCompareColor,
  ]);

  // Recursively try to find the parent row element with the
  // group data attribute
  const MAX_RECURSION_DEPTH = 10;
  const getRowGroup = (node, level = 1) => {
    // If we end up hitting the root of the DOM tree,
    // or if we exhaust the recursion depth,
    // then return undefined
    if (node === null || level > MAX_RECURSION_DEPTH) {
      return undefined;
    }

    const group = node.getAttribute('data-group');
    if (group === null) {
      return getRowGroup(node.parentElement, (level = level + 1));
    } else {
      return group;
    }
  };

  const onTableHover = (e) => {
    // console.log(e);
    const group = getRowGroup(e.target);

    // We didn't find a parent row group,
    // probably because we're hovering over the header
    // or some other table element
    if (group === undefined) {
      updateHoverGroup(null);
    }

    updateHoverGroup(group);
  };

  const updateHoverGroup = _.debounce((hoverGroup) => {
    // Don't fire the action if there's no change
    if (hoverGroup === configStore.hoverGroup) {
      return;
    }
    configStore.updateHoverGroup(hoverGroup);
  }, 20);

  const onRowClick = (rowIndex, row) => {
    //console.log(rowIndex, row, column);

    let newGroups;

    let selectedGroup = row.group;

    // If we selected the reference group, and we're in lineage/clade mode,
    // then ignore
    if (
      selectedGroup === GROUPS.REFERENCE_GROUP &&
      Object.keys(config.group_cols).includes(configStore.groupKey)
    ) {
      return;
    }

    // If the item is already selected, then deselect it
    if (
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

  const onSort = (sortColumn, sortDirection) => {
    sortDirection =
      sortDirection === SORT_DIRECTIONS.SORT_NONE
        ? SORT_DIRECTIONS.SORT_ASC
        : sortDirection;
    // console.log(sortColumn, sortDirection);

    plotSettingsStore.setTableSort(sortColumn, sortDirection);
  };

  useEffect(() => {
    const rows = sortRows(
      state.rows,
      comparer({
        sortDirection: plotSettingsStore.tableSortDirection,
        sortColumn: plotSettingsStore.tableSortColumn,
      })
    );

    setState({
      ...state,
      rows,
    });
  }, [plotSettingsStore.tableSortColumn, plotSettingsStore.tableSortDirection]);

  if (UIStore.caseDataState === ASYNC_STATES.STARTED) {
    return (
      <div
        style={{
          paddingRight: '24px',
          paddingLeft: '12px',
          paddingTop: '24px',
          height: '100%',
        }}
      >
        {_.times(20, (i) => (
          <SkeletonElement
            key={Math.random()}
            delay={5 + i + (i % 2) * 12.5}
            height={25}
          />
        ))}
      </div>
    );
  }

  // If we have no rows, then return an empty element
  // We'll always have the "reference" row, so no rows = 1 row
  if (
    state.rows.length === 1 ||
    configStore.selectedLocationNodes.length === 0
  ) {
    return <EmptyDataTable />;
  }

  //console.log(state.rows);
  return (
    <DataTableContainer>
      <TableOptions />
      <span className="position-title" style={{ marginLeft: 220 }}>
        {configStore.dnaOrAa === DNA_OR_AA.DNA
          ? 'Genomic Coordinate'
          : 'Residue Index'}
      </span>
      <div style={{ paddingLeft: '10px' }} onMouseMove={onTableHover}>
        <DataTable
          posColOffset={4}
          columns={state.columns}
          rows={state.rows}
          rowGetter={(i) => state.rows[i]}
          rowsCount={state.rows ? state.rows.length : 0}
          height={state.tableHeight}
          headerRowHeight={45}
          filterRowHeight={45}
          rowHeight={25}
          minColumnWidth={25}
          sortColumn={plotSettingsStore.tableSortColumn}
          sortDirection={plotSettingsStore.tableSortDirection}
          onSort={onSort}
          rowRenderer={RowRenderer}
          onRowClick={onRowClick}
          enableCellSelect={false}
          enableCellAutoFocus={false}
        />
      </div>
    </DataTableContainer>
  );
});

export default NewLineageDataTable;
