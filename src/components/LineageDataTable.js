/* eslint-disable react/display-name */
import React from 'react';

import _ from 'underscore';

import DataTable from 'react-data-table-component';
import HeatmapCell from './Cells/HeatmapCell';
import AddToSidepanelCheckbox from './AddToSidepanelCheckbox';
import { useStores } from '../stores/connect';
import { getReferenceSequence } from '../utils/lineageData';

const LineageDataTable = () => {
  const { covidStore } = useStores();

  // Get the bases at the positions, for the reference sequence
  let ref_seq = getReferenceSequence();

  // Build a column for each changing position
  let pos_cols = [];
  covidStore.changingPositions.forEach((pos) => {
    pos_cols.push({
      name: pos.toString(),
      selector: 'pos_' + pos.toString(),
      sortable: false,
      width: '40px',
      center: true,
      compact: true,
      style: {
        fontFamily: 'monospace',
        fontWeight: '500',
        fontSize: '1.25em',
      },
      conditionalCellStyles: [
        {
          when: (row) => row['pos_' + pos.toString()] != ref_seq[pos],
          style: {
            backgroundColor: '#FFFF00',
          },
        },
      ],
    });
  });

  let maxCasesSum = _.reduce(
    covidStore.caseDataAggLineageList,
    (memo, lineage) => Math.max(memo, lineage.cases_sum),
    0
  );
  let minCasesSum = _.reduce(
    covidStore.caseDataAggLineageList,
    (memo, lineage) => Math.min(memo, lineage.cases_sum),
    0
  );
  let maxCasesPercent = _.reduce(
    covidStore.caseDataAggLineageList,
    (memo, lineage) => Math.max(memo, lineage.cases_percent),
    0
  );
  let minCasesPercent = _.reduce(
    covidStore.caseDataAggLineageList,
    (memo, lineage) => Math.min(memo, lineage.cases_percent),
    0
  );

  console.log(covidStore);

  return (
    <DataTable
      className="data-table"
      data={covidStore.caseDataAggLineageList}
      columns={[
        {
          name: 'Lineage',
          selector: 'lineage',
          sortable: true,
          width: '100px',
          style: {
            fontWeight: '700',
          },
        },
        {
          name: 'Cases',
          selector: 'cases_sum',
          sortable: true,
          width: '85px',
          cell: (row) => {
            return (
              <HeatmapCell
                value={row.cases_sum}
                min={minCasesSum}
                max={maxCasesSum}
                percent={false}
              />
            );
          },
        },
        {
          selector: 'cases_percent',
          sortable: true,
          width: '85px',
          cell: (row) => {
            return (
              <HeatmapCell
                value={row.cases_percent}
                min={minCasesPercent}
                max={maxCasesPercent}
                percent={true}
              />
            );
          },
        },
        {
          name: 'is in sidepanel',
          selector: null,
          sortable: false,
          width: '100%',
          cell: () => {
            return <AddToSidepanelCheckbox />;
          },
        },
      ].concat(pos_cols)}
      striped={true}
      highlightOnHover={true}
      dense={true}
      // fixedHeader={true}
      // fixedHeaderScrollHeight={'400px'}

      pagination={false}
      defaultSortField={'lineage'}
      defaultSortAsc={true}
      conditionalRowStyles={[
        {
          when: (row) => row.lineage == 'Reference',
          style: 'background-color: #dff3fe !important;',
        },
      ]}
      sortFunction={(rows, field, direction) => {
        // Set aside the reference, and remove it from the rows list
        let refRow = _.findWhere(rows, { lineage: 'Reference' });
        rows = _.reject(rows, (row) => row.lineage == 'Reference');

        // Normal sorting...
        rows = _.sortBy(rows, (row) => {
          return row[field];
        });
        // Reverse if descending
        if (direction == 'desc') {
          rows.reverse();
        }
        // Add the reference row to the beginning
        rows.unshift(refRow);

        return rows;
      }}
      customStyles={{
        headCells: {
          style: {
            paddingLeft: '8px', // override the cell padding for head cells
            paddingRight: '8px',
          },
        },
        cells: {
          style: {
            paddingLeft: '8px', // override the cell padding for data cells
            paddingRight: '8px',
          },
        },
      }}
    />
  );
};

LineageDataTable.displayName = 'LineageDataTable';

export default LineageDataTable;
