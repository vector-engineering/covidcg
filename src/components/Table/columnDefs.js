/* eslint-disable react/display-name */
import React from 'react';
import LetterCell from '../Cells/LetterCell';
import HeatmapCell from '../Cells/HeatmapCell';
import PosHeaderCell from '../Cells/PosHeaderCell';
import AddToSidepanelCheckbox from '../AddToSidepanelCheckbox';
import { snapGeneHighlightColors } from '../../utils/colors';

export const positionColumn = () => ({
  name: 'Position',
  key: 'pos',
  sortable: true,
  width: 70,
  cellClass: 'no-overflow',
  frozen: true,
});

export const indexColumn = () => ({
  name: 'Index',
  key: 'pos',
  sortable: true,
  width: 50,
  frozen: true,
});
export const refColumn = () => ({
  name: 'Ref',
  key: 'ref',
  sortable: true,
  width: 36,
  frozen: true,
  cellClass: 'no-padding',
  formatter: (val) => <LetterCell value={val.row['ref']} />,
});
export const altColumn = () => ({
  name: 'Alt',
  key: 'alt',
  sortable: true,
  width: 36,
  frozen: true,
  cellClass: 'no-padding',
  formatter: (val) => <LetterCell value={val.row['alt']} />,
});

export const geneColumn = () => ({
  name: 'Gene',
  key: 'gene',
  sortable: true,
  width: 65,
  cellClass: 'no-overflow',
  frozen: true,
});

export const lineageColumn = () => ({
  key: 'group',
  name: 'Lineage',
  selector: 'group',
  sortable: true,
  width: 85,
  frozen: true,
});

const conditionCompare = (base, refBase, matchOrMismatch) => {
  // Flip the XOR (XNOR)
  return !((base === refBase) ^ (matchOrMismatch === 'match' ? true : false));
};

export const getSinglePosColumn = ({
  pos,
  col,
  colorMode,
  refRow,
  compareMode,
  compareColor,
  colors,
}) => ({
  name: pos.toString(),
  key: col,
  cellClass: 'pos',
  sortable: false,
  width: 25,
  formatter: (val) => {
    // console.log(val);
    const row = val.row;
    let cellBgColor = 'transparent';
    // Define the coloring behavior
    if (
      colorMode === 'compare' &&
      conditionCompare(row[col], refRow[col], compareMode)
    ) {
      cellBgColor =
        compareColor === 'code'
          ? colors[row[col]]
          : snapGeneHighlightColors[compareColor];
    }
    // Just coloring by code
    else if (colorMode === 'code') {
      cellBgColor = colors[row[col]];
    }

    return <LetterCell value={row[col]} bgColor={cellBgColor} />;
  },
  headerRenderer: (val) => {
    return <PosHeaderCell pos={val.column.name} />;
  },
});

export const getDefaultColumns = ({
  minCasesPercent,
  maxCasesPercent,
  minCasesSum,
  maxCasesSum,
}) => [
  {
    name: 'Seqs',
    key: 'cases_sum',
    sortable: true,
    width: 55,
    frozen: true,
    cellClass: 'no-padding',
    formatter: (val) => {
      const row = val.row;
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
    name: '% Seqs',
    key: 'cases_percent',
    sortable: true,
    width: 70,
    frozen: true,
    cellClass: 'no-padding',
    formatter: (val) => {
      const row = val.row;
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
    name: 'Show structure',
    key: null,
    sortable: false,
    width: 40,
    frozen: true,
    formatter: (val) => <AddToSidepanelCheckbox groupKey={val.row.group} />,
  },
];
