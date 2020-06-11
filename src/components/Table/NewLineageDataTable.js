/* eslint-disable react/display-name */
import React, { useState, useEffect } from 'react';
import styled from 'styled-components';

import DataGrid from 'react-data-grid';
import 'react-data-grid/dist/react-data-grid.css';
import { observer } from 'mobx-react';
import _ from 'underscore';

import { useStores } from '../../stores/connect';
import { nanmax, nanmin } from '../../utils/math';
import { snapGeneNTColors, shingAAColors } from '../../utils/colors';
import TableOptions from './TableOptions';
import {
  geneColumn,
  positionColumn,
  indexColumn,
  refColumn,
  altColumn,
  lineageColumn,
  getDefaultColumns,
  getSinglePosColumn,
} from './columnDefs';

const DataTableContainer = styled.div`
  span.position-title {
    font-size: 12px;
    font-weight: 500;
  }
`;

const comparer = ({ sortDirection, sortColumn }) => (a, b) => {
  if (sortDirection === 'ASC') {
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
  const { covidStore } = useStores();

  // Color by 'compare': Comparison to reference, or 'code': With a defined color code
  const [state, setState] = useState({
    colorMode: 'compare',
    compareMode: 'mismatch',
    compareColor: 'yellow',
    rows: covidStore.caseDataAggGroup,
    sortColumn: 'group',
    sortDirection: 'ASC',
  });

  useEffect(() => {
    setState({
      ...state,
      rows: sortRows(
        covidStore.caseDataAggGroup,
        comparer({
          sortDirection: state.sortDirection,
          sortColumn: state.sortColumn,
        })
      ),
    });
  }, [covidStore.caseDataAggGroup]);

  const handleColorModeChange = (event) =>
    setState({ ...state, colorMode: event.target.value });
  const handleCompareModeChange = (event) =>
    setState({ ...state, compareMode: event.target.value });
  const handleCompareColorChange = (event) =>
    setState({ ...state, compareColor: event.target.value });

  // Get the maximum and minimum cases_sum and cases_percent for the colormaps
  // Ignore those values for the reference row (which are NaN)
  let maxCasesSum = _.reduce(
    covidStore.caseDataAggGroup,
    (memo, group) => nanmax(memo, group.cases_sum),
    0
  );
  let minCasesSum = _.reduce(
    covidStore.caseDataAggGroup,
    (memo, group) => nanmin(memo, group.cases_sum),
    0
  );
  let maxCasesPercent = _.reduce(
    covidStore.caseDataAggGroup,
    (memo, group) => nanmax(memo, group.cases_percent),
    0
  );
  let minCasesPercent = _.reduce(
    covidStore.caseDataAggGroup,
    (memo, group) => nanmin(memo, group.cases_percent),
    0
  );

  const handleGridSort = (sortColumn, sortDirection) => {
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
    if (covidStore.groupKey === 'lineage') {
      _columns.push(lineageColumn(handleGridSort));
    }

    // For SNP grouping, add each SNP chunk as its own column
    if (covidStore.groupKey === 'snp') {
      // Add the gene column, if we're in AA mode
      if (covidStore.dnaOrAa === 'aa') {
        _columns.push(geneColumn(handleGridSort));
      }
      // Add the position column
      // We don't need as much space for this, for AA mode
      if (covidStore.dnaOrAa === 'dna') {
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
    let refRow = _.findWhere(covidStore.caseDataAggGroup, {
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

      let colors =
        covidStore.dnaOrAa === 'dna' ? snapGeneNTColors : shingAAColors;

      // 0-indexed to 1-indexed
      let pos = parseInt(col.substring(4)) + 1;

      _columns.push(
        getSinglePosColumn({
          pos,
          col,
          colorMode: state.colorMode,
          refRow,
          compareMode: state.compareColor,
          compareColor: state.compareColor,
          colors,
        })
      );
    });
    return _columns;
  };

  const columns = buildColumns() || [];

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
      <span className="position-title">
        {covidStore.dnaOrAa === 'dna' ? 'Genomic Coordinate' : 'Residue Index'}
      </span>
      <DataGrid
        columns={columns}
        rowGetter={(i) => state.rows[i]}
        rows={state.rows}
        rowsCount={state.rows.length}
        minHeight={500}
        sortColumn={state.sortColumn}
        sortDirection={state.sortDirection}
        onSort={handleGridSort}
      />
    </DataTableContainer>
  );
});

export default NewLineageDataTable;
