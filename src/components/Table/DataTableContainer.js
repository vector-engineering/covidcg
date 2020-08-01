/* eslint-disable react/display-name */
import React, { useState, useEffect } from 'react';
import styled from 'styled-components';

import EmptyDataTable from './EmptyDataTable';
import 'react-data-grid/dist/react-data-grid.css';
import { observer } from 'mobx-react';
import _ from 'underscore';

import { useStores } from '../../stores/connect';
import { nanmax, nanmin } from '../../utils/math';
import {
  snapGeneNTColors,
  shingAAColors,
  clustalXAAColors,
  zappoAAColors,
  transmembraneAAColors,
} from '../../constants/colors';
import TableOptions from './TableOptions';
import {
  geneColumn,
  proteinColumn,
  positionColumn,
  indexColumn,
  refColumn,
  altColumn,
  lineageColumn,
  cladeColumn,
  getDefaultColumns,
  getSinglePosColumn,
} from './columnDefs';
import SkeletonElement from '../Common/SkeletonElement';
import { asyncStates } from '../../stores/UIStore';
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
  let refRow = _.findWhere(rows, { group: 'Reference' });
  rows = _.reject(rows, (row) => row.group == 'Reference');
  rows = rows.sort(sortFn);
  rows.unshift(refRow);
  return rows;
};

const NewLineageDataTable = observer(() => {
  const { dataStore, UIStore, configStore } = useStores();

  const [state, setState] = useState({
    // Color by 'compare': Comparison to reference, or 'code': With a defined color code
    colorMode: 'compare',
    // 'match' or 'mismatch'
    compareMode: 'mismatch',
    compareColor: 'yellow',
    rows: dataStore.caseDataAggGroup,
    sortColumn: 'cases_sum',
    sortDirection: 'DESC',
  });

  // Reset the coloring settings when dnaOrAa is changed
  useEffect(() => {
    setState({
      ...state,
      colorMode: 'compare',
      compareMode: 'mismatch',
      compareColor: 'yellow',
    });
  }, [configStore.dnaOrAa]);

  useEffect(() => {
    setState({
      ...state,
      rows: sortRows(
        dataStore.caseDataAggGroup,
        comparer({
          sortDirection: state.sortDirection,
          sortColumn: state.sortColumn,
        })
      ),
    });
  }, [dataStore.caseDataAggGroup]);

  const handleColorModeChange = (event) =>
    setState({ ...state, colorMode: event.target.value });
  const handleCompareModeChange = (event) =>
    setState({ ...state, compareMode: event.target.value });
  const handleCompareColorChange = (event) =>
    setState({ ...state, compareColor: event.target.value });

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

  const updateHoverGroup = (hoverGroup) => {
    // Don't fire the action if there's no change
    if (hoverGroup === configStore.hoverGroup) {
      return;
    }
    configStore.updateHoverGroup(hoverGroup);
  };

  const onRowClick = (rowIndex, row) => {
    //console.log(rowIndex, row, column);

    let newGroups;

    let selectedGroup = row.group;
    // If we selected a group that's grouped into 'other',
    // then pretend like we're selecting 'other' instead
    if (!dataStore.groupsToKeep.includes(row.group)) {
      selectedGroup = 'other';
    }

    // If we selected the reference group, and we're in lineage/clade mode,
    // then ignore
    if (
      selectedGroup === 'Reference' &&
      (configStore.groupKey === 'lineage' || configStore.groupKey === 'clade')
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

  const renderTable = () => {
    if (
      UIStore.caseDataState === asyncStates.STARTED ||
      UIStore.aggCaseDataState === asyncStates.STARTED
    ) {
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
    if (state.rows.length === 1) {
      return <EmptyDataTable />;
    }

    // Get the maximum and minimum cases_sum and cases_percent for the colormaps
    // Ignore those values for the reference row (which are NaN)
    let maxCasesSum = _.reduce(
      dataStore.caseDataAggGroup,
      (memo, group) => nanmax(memo, group.cases_sum),
      0
    );
    let minCasesSum = _.reduce(
      dataStore.caseDataAggGroup,
      (memo, group) => nanmin(memo, group.cases_sum),
      0
    );
    let maxCasesPercent = _.reduce(
      dataStore.caseDataAggGroup,
      (memo, group) => nanmax(memo, group.cases_percent),
      0
    );
    let minCasesPercent = _.reduce(
      dataStore.caseDataAggGroup,
      (memo, group) => nanmin(memo, group.cases_percent),
      0
    );

    const handleGridSort = (sortColumn, sortDirection) => {
      // console.log('handle grid sort', sortColumn, sortDirection);

      let _sortDirection = sortDirection;
      if (sortDirection === 'NONE') {
        _sortDirection = 'ASC';
      }

      const _rows = sortRows(
        state.rows,
        comparer({
          sortDirection: _sortDirection,
          sortColumn: state.sortColumn,
        })
      );

      setState({
        ...state,
        rows: _rows,
        sortColumn,
        sortDirection: _sortDirection,
      });
    };

    const buildColumns = () => {
      let _columns = [];
      // For lineage grouping, add lineage column
      if (configStore.groupKey === 'lineage') {
        _columns.push(lineageColumn(handleGridSort));
      } else if (configStore.groupKey === 'clade') {
        _columns.push(cladeColumn(handleGridSort));
      }

      // For SNP grouping, add each SNP chunk as its own column
      if (configStore.groupKey === 'snp') {
        // Add the gene column, if we're in AA mode
        if (configStore.dnaOrAa === 'aa') {
          if (configStore.coordinateMode === 'gene') {
            _columns.push(geneColumn(handleGridSort));
          } else if (configStore.coordinateMode === 'protein') {
            _columns.push(proteinColumn(handleGridSort));
          }
        }
        // Add the position column
        // We don't need as much space for this, for AA mode
        if (configStore.dnaOrAa === 'dna') {
          _columns.push(positionColumn(handleGridSort));
        } else {
          _columns.push(indexColumn(handleGridSort));
        }
        _columns.push(refColumn(handleGridSort));
        _columns.push(altColumn(handleGridSort));
      }

      // example row
      // cases_percent: 0.6880290205562273
      // cases_sum: 569
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
          handleGridSort,
        })
      );

      // Build a column for each changing position
      let refRow = _.findWhere(dataStore.caseDataAggGroup, {
        group: 'Reference',
      });
      if (!refRow) {
        return null;
      }

      Object.keys(refRow).forEach((col) => {
        // Only process columns starting with "pos_"
        if (!col.startsWith('pos_')) {
          return;
        }

        let colors;
        if (configStore.dnaOrAa === 'dna') {
          colors = snapGeneNTColors;
        } else {
          if (state.compareColor === 'code' || state.colorMode === 'code') {
            colors = shingAAColors;
          } else if (
            state.compareColor === 'clustal' ||
            state.colorMode === 'clustal'
          ) {
            colors = clustalXAAColors;
          } else if (
            state.compareColor === 'zappo' ||
            state.colorMode === 'zappo'
          ) {
            colors = zappoAAColors;
          } else if (
            state.compareColor === 'zhao-london' ||
            state.colorMode === 'zhao-london'
          ) {
            colors = transmembraneAAColors;
          }
        }

        // 0-indexed to 1-indexed
        let pos = parseInt(col.substring(4));
        if (configStore.dnaOrAa === 'dna') {
          pos += 1;
        }
        if (configStore.groupKey === 'snp' && configStore.dnaOrAa === 'aa') {
          pos += 1;
        }

        _columns.push(
          getSinglePosColumn({
            pos,
            col,
            colorMode: state.colorMode,
            refRow,
            compareMode: state.compareMode,
            compareColor: state.compareColor,
            colors,
          })
        );
      });

      // console.log(_columns);

      return _columns;
    };

    const columns = buildColumns() || [];

    let positionTitleOffset = 0;
    let posColOffset = 0;
    if (
      configStore.groupKey === 'lineage' ||
      configStore.groupKey === 'clade'
    ) {
      positionTitleOffset = 220;
      posColOffset = 4;
    } else if (configStore.groupKey === 'snp') {
      if (configStore.dnaOrAa === 'dna') {
        positionTitleOffset = 280;
        posColOffset = 6;
      } else {
        positionTitleOffset = 325;
        posColOffset = 7;
      }
    }

    // console.log(state.sortColumn, state.sortDirection, columns);
    // console.log(state.rows);

    return (
      <>
        <span
          className="position-title"
          style={{ marginLeft: positionTitleOffset }}
        >
          {configStore.dnaOrAa === 'dna'
            ? 'Genomic Coordinate'
            : 'Residue Index'}
        </span>
        <div style={{ paddingLeft: '10px' }} onMouseMove={onTableHover}>
          <DataTable
            posColOffset={posColOffset}
            columns={columns}
            rows={state.rows}
            rowGetter={(i) => state.rows[i]}
            rowsCount={state.rows ? state.rows.length : 0}
            height={state.tableHeight}
            headerRowHeight={45}
            filterRowHeight={45}
            rowHeight={25}
            minColumnWidth={25}
            sortColumn={state.sortColumn}
            sortDirection={state.sortDirection}
            onSort={handleGridSort}
            rowRenderer={RowRenderer}
            onRowClick={onRowClick}
          />
        </div>
      </>
    );
  };

  return (
    <DataTableContainer>
      <TableOptions
        handleColorModeChange={handleColorModeChange}
        handleCompareModeChange={handleCompareModeChange}
        handleCompareColorChange={handleCompareColorChange}
        colorMode={state.colorMode}
        compareColor={state.compareColor}
        compareMode={state.compareMode}
      />
      {renderTable()}
    </DataTableContainer>
  );
});

export default NewLineageDataTable;
