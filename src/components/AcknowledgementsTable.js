import React, { useState, useEffect } from 'react';
import styled from 'styled-components';
import { observer } from 'mobx-react';
import _ from 'underscore';

import { useStores } from '../stores/connect';

import DataGrid from 'react-data-grid';
import AckCell from './Cells/AckCell';
import AckAuthorCell from './Cells/AckAuthorCell';

import { getAckTextsFromAckIds } from '../utils/acknowledgements';

const AckContainer = styled.div`
  /* Data grid styles */
  .rdg {
    font-size: 12px;
  }

  /* All cells */
  .rdg-cell {
    // Allow text selection
    user-select: text;
  }

  .rdg-header-row {
  }

  /* All other rows */
  .rdg-row {
    font-weight: normal;
  }
`;

const comparer = ({ sortDirection, sortColumn }) => (a, b) => {
  if (sortDirection === 'ASC') {
    return a[sortColumn] > b[sortColumn] ? 1 : -1;
  } else if (sortDirection === 'DESC') {
    return a[sortColumn] < b[sortColumn] ? 1 : -1;
  }
};

const sortRows = (rows, sortFn) => {
  return rows.sort(sortFn);
};

const AcknowledgementsTable = observer(() => {
  const { covidStore } = useStores();

  const [state, setState] = useState({
    rows: [],
    sortColumn: 'cases_sum',
    sortDirection: 'DESC',
  });

  useEffect(() => {
    let ackIds = _.pluck(covidStore.selectedRows, 'ack_id');
    // Get unique acknowledgement IDs
    ackIds = Array.from(new Set(ackIds));
    // Remove null IDs
    ackIds = _.reject(ackIds, (ackId) => ackId == -1);

    // Get the list of selected Accession IDs, and map to
    // acknowledgement texts
    let ackTexts = getAckTextsFromAckIds(ackIds);
    // Set the acknowledgement ID for each acknowledgement object
    for (let i = 0; i < ackTexts.length; i++) {
      ackTexts[i]['ack_id'] = ackIds[i];
    }

    setState({
      ...state,
      rows: sortRows(
        ackTexts,
        comparer({
          sortDirection: state.sortDirection,
          sortColumn: state.sortColumn,
        })
      ),
    });
  }, [covidStore.caseDataAggGroup]);

  const handleGridSort = (sortColumn, sortDirection) => {
    let _sortDirection = sortDirection;

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

  console.log(state.rows);

  return (
    <AckContainer>
      <DataGrid
        columns={[
          {
            name: 'Originating Lab',
            key: 'Originating lab',
            sortable: true,
            /* eslint-disable react/display-name */
            formatter: (val) => {
              return <AckCell text={val.row[val.column.key]} />;
            },
          },
          {
            name: 'Submitting Lab',
            key: 'Submitting lab',
            sortable: true,
            /* eslint-disable react/display-name */
            formatter: (val) => {
              return <AckCell text={val.row[val.column.key]} />;
            },
          },
          {
            name: 'Authors',
            key: 'authors_abbrev',
            sortable: true,
            /* eslint-disable react/display-name */
            formatter: (val) => {
              return (
                <AckAuthorCell
                  shortText={val.row['authors_abbrev']}
                  longText={val.row['Authors']}
                />
              );
            },
          },
        ]}
        rows={state.rows}
        rowsCount={state.rows.length}
        rowGetter={(i) => state.rows[i]}
        rowKey="ack_id"
        height={300}
        headerRowHeight={30}
        filterRowHeight={30}
        rowHeight={36}
        minColumnWidth={25}
        sortColumn={state.sortColumn}
        sortDirection={state.sortDirection}
        onSort={handleGridSort}
      />
    </AckContainer>
  );
});

export default AcknowledgementsTable;
