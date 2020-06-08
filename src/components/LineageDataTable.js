/* eslint-disable react/display-name */
import React from 'react';
import styled from 'styled-components';
import _ from 'underscore';

import DataTable from 'react-data-table-component';
import HeatmapCell from './Cells/HeatmapCell';
import AddToSidepanelCheckbox from './AddToSidepanelCheckbox';
import { useStores } from '../stores/connect';
import { nanmin, nanmax } from '../utils/math';

const StyledDataTable = styled(DataTable)`
  // hide the header
  header {
    display: none;
  }

  .rdt_TableHeadRow {
    padding-top: 10px;

    .rdt_TableCol {
      padding-left: 4px;
      padding-right: 4px;
    }
  }

  .rdt_TableCol_Sortable[id^='column-pos'] {
    // color: red;
    & > div {
      margin-left: 15px;
      margin-bottom: 10px;
      transform: rotate(-45deg);
    }
  }

  .rdt_TableRow {
    .rdt_TableCell {
      padding-left: 4px;
      padding-right: 4px;
    }
  }
`;

const sortFn = (rows, field, direction) => {
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
};

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
      width: '24px',
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

  // Better column names
  let groupKeyName = '';
  if (covidStore.groupKey === 'lineage') {
    groupKeyName = 'Lineage';
  } else if (covidStore.groupKey === 'snp') {
    if (covidStore.dnaOrAa === 'dna') {
      groupKeyName = 'NT SNP';
    } else {
      groupKeyName = 'AA SNP';
    }
  }

  return (
    <StyledDataTable
      className="data-table"
      data={covidStore.caseDataAggGroup}
      keyField="group"
      columns={[
        {
          name: groupKeyName,
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
          width: '60px',
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
          width: '75px',
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
          width: '80px',
          cell: (row) => {
            return <AddToSidepanelCheckbox groupKey={Math.random()} />;
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
      sortFunction={sortFn}
    />
  );
};

LineageDataTable.displayName = 'LineageDataTable';

export default LineageDataTable;
