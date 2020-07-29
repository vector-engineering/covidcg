import React, { useCallback, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import DataGrid from 'react-data-grid';

const DataGridContainer = styled.div`
  margin-bottom: 24px;
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
    // background-color: transparent;
    line-height: 65px;
    height: 45px;
    font-weight: 500;
    z-index: 4;
    width: 100%;

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
      .rdg-header-sort-cell {
        display: block;
        position: relative;
        height: 100%;
        // Sorting caret
        span:nth-child(2) {
          position: absolute;
          line-height: normal;
          top: 5px;
          left: 4px;
        }
      }
    }
  }

  /* All other rows */
  .rdg-row {
    background-color: transparent;
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

const DataTable = ({ posColOffset, rows, ...rest }) => {
  const [ref] = useHookWithRefCallback();

  return (
    <DataGridContainer posColOffset={posColOffset}>
      <div ref={ref} />
      <DataGrid rows={rows} height={400} {...rest} />
    </DataGridContainer>
  );
};

DataTable.propTypes = {
  posColOffset: PropTypes.number,
  rows: PropTypes.array.isRequired,
};

export default DataTable;
