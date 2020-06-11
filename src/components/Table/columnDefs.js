/* eslint-disable react/display-name */
import React from 'react';
import LetterCell from '../Cells/LetterCell';
import HeatmapCell from '../Cells/HeatmapCell';
import AddToSidepanelCheckbox from '../AddToSidepanelCheckbox';
import { snapGeneHighlightColors } from '../../utils/colors';

export const positionColumn = () => ({
  name: 'Position',
  key: 'pos',
  sortable: true,
  width: 60,
});

export const indexColumn = () => ({
  name: 'Index',
  key: 'pos',
  sortable: true,
  width: 45,
});
export const refColumn = () => ({
  name: 'Ref',
  key: 'ref',
  sortable: true,
  width: 24,
  formatter: (val) => <LetterCell value={val.row['ref']} />,
});
export const altColumn = () => ({
  name: 'Alt',
  key: 'alt',
  sortable: true,
  width: 24,
  formatter: (val) => <LetterCell value={val.row['alt']} />,
});

export const geneColumn = () => ({
  name: 'Gene',
  key: 'gene',
  sortable: true,
  width: 50,
});

export const lineageColumn = () => ({
  key: 'group',
  name: 'Lineage',
  selector: 'group',
  sortable: true,
  width: 100,
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
  sortable: false,
  formatter: (val) => {
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
    width: 60,
    formatter: (val) => {
      const row = val.row;
      return (
        <HeatmapCell
          value={row.cases_sum}
          min={minCasesSum}
          max={maxCasesSum}
          percent={true}
        />
      );
    },
  },
  {
    name: '% Seqs',
    key: 'cases_percent',
    sortable: true,
    width: 75,
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
    width: 80,
    formatter: (val) => <AddToSidepanelCheckbox groupKey={val.row.group} />,
  },
];
