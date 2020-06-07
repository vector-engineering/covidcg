/* eslint-disable react/display-name */
import React from 'react';

import _ from 'underscore';

import DataTable from 'react-data-table-component';
import HeatmapCell from './Cells/HeatmapCell';
import AddToSidepanelCheckbox from './AddToSidepanelCheckbox';
import { useStores } from '../stores/connect';
import { nanmin, nanmax } from '../utils/math';

const LineageDataTable = () => {
  const { covidStore } = useStores();

  // Build a column for each changing position
  let posCols = [];
  let refRow = _.findWhere(covidStore.caseDataAggGroup, {
    group: 'Reference',
  });
  Object.keys(refRow).forEach((col) => {
    // Only process columns starting with "pos_"
    if (!col.startsWith('pos_')) {
      return;
    }

    // 0-indexed to 1-indexed
    let pos = parseInt(col.substring(4)) + 1;

    posCols.push({
      name: pos.toString(),
      selector: col,
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
          when: (row) => row[col] != refRow[col],
          style: {
            backgroundColor: '#FFFF00',
          },
        },
      ],
    });
  });

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

  console.log(maxCasesPercent, minCasesPercent);

  // console.log(covidStore.caseDataAggGroup)

  return (
    <DataTable
      className="data-table"
      data={covidStore.caseDataAggGroup}
      columns={[
        {
          name: 'Lineage',
          selector: 'group',
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
          name: '% Cases',
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
      ].concat(posCols)}
      striped={true}
      highlightOnHover={true}
      dense={true}
      // fixedHeader={true}
      // fixedHeaderScrollHeight={'400px'}

      pagination={false}
      defaultSortField={'group'}
      defaultSortAsc={true}
      conditionalRowStyles={[
        {
          when: (row) => row.group == 'Reference',
          style: 'background-color: #dff3fe !important;',
        },
      ]}
      sortFunction={(rows, field, direction) => {
        // Set aside the reference, and remove it from the rows list
        let refRow = _.findWhere(rows, { group: 'Reference' });
        rows = _.reject(rows, (row) => row.group == 'Reference');

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
