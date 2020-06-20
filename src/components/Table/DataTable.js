import React, { useLayoutEffect, useCallback, useState, useRef } from 'react';
import styled from 'styled-components';
import DataGrid from 'react-data-grid';
import { useStores } from '../../stores/connect';

const DataGridContainer = styled.div`
  /* Data grid styles */
  .rdg {
    border-top: none;
    border-left: none;
    font-size: 12px;
  }

  /* All cells */
  .rdg-cell {
    // Allow text selection
    user-select: text;

    &.no-padding {
      padding: 0px;
    }
    &.no-overflow {
      text-overflow: unset;
    }
  }

  .rdg-header-row {
    background-color: white;
    line-height: 65px;
    height: 45px;
    font-weight: 500;

    .rdg-cell {
      font-size: 12px;
      border-right: none;
      &.rdg-cell-frozen-last {
        box-shadow: none;
      }
    }
    /* Position columns */
    .rdg-cell:nth-child(n + ${(props) => props.posColOffset}) {
      padding: 0px;
      background-color: transparent;
      overflow: unset;
    }
  }

  /* All other rows */
  .rdg-row {
    .rdg-cell {
      border-right: none;
      &.rdg-cell-frozen-last {
        border-right: 1px solid #ccc;
        box-shadow: none;
      }
    }
    /* Position cells */
    .rdg-cell.pos {
      padding: 0px;
    }
  }
`;

DataGridContainer.defaultProps = {
  posColOffset: 0,
};

function useHookWithRefCallback() {
  const ref = useRef(null);
  const setRef = useCallback((node) => {
    if (ref.current) {
      // Make sure to cleanup any events/references added to the last instance
    }

    if (node) {
      // Check if a node is actually passed. Otherwise node would be null.
      // You can now do what you need to, addEventListeners, measure, etc.
    }

    // Save a reference to the node
    ref.current = node;
  }, []);

  return [setRef];
}

const DataTable = ({
  posColOffset,
  rows,
  columns,
  sortColumn,
  sortDirection,
  handleGridSort,
}) => {
  const [ref] = useHookWithRefCallback();

  return (
    <DataGridContainer posColOffset={posColOffset}>
      <div ref={ref} />
      <DataGrid
        columns={columns}
        rowGetter={(i) => rows[i]}
        rows={rows}
        rowsCount={rows ? rows.length : 0}
        height={400}
        headerRowHeight={45}
        filterRowHeight={45}
        rowHeight={25}
        minColumnWidth={25}
        sortColumn={sortColumn}
        sortDirection={sortDirection}
        onSort={handleGridSort}
      />
    </DataGridContainer>
  );
};

export default DataTable;
