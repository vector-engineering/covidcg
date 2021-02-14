/* eslint-disable react/display-name */
import React from 'react';
import GroupCell from '../Cells/GroupCell';
import LetterCell from '../Cells/LetterCell';
import HeatmapCell from '../Cells/HeatmapCell';
import PosHeaderCell from '../Cells/PosHeaderCell';
// import AddToSidepanelCheckbox from '../AddToSidepanelCheckbox';
import { snapGeneHighlightColors } from '../../constants/colors';
import {
  GROUPS,
  COLOR_MODES,
  COMPARE_MODES,
  COMPARE_COLORS,
} from '../../constants/defs.json';

export const positionColumn = () => ({
  name: 'Position',
  key: 'pos',
  sortable: true,
  width: 70,
  cellClass: 'no-overflow',
  frozen: true,
  formatter: (val) => (
    <GroupCell text={val.row['pos']} color={val.row['color']} />
  ),
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
  width: 120,
  cellClass: 'no-overflow',
  frozen: true,
  formatter: (val) => (
    <GroupCell text={val.row['gene']} color={val.row['color']} />
  ),
});

export const proteinColumn = () => ({
  name: 'Protein',
  key: 'protein',
  sortable: true,
  width: 120,
  cellClass: 'no-overflow',
  frozen: true,
  formatter: (val) => (
    <GroupCell text={val.row['protein']} color={val.row['color']} />
  ),
});

export const groupColumn = ({ title }) => ({
  key: 'group',
  name: title,
  selector: 'group',
  sortable: true,
  width: 85,
  frozen: true,
  formatter: (val) => (
    <GroupCell text={val.row['group']} color={val.row['color']} />
  ),
});

const conditionCompare = (base, refBase, matchOrMismatch) => {
  // Flip the XOR (XNOR)
  return !(
    (base === refBase) ^
    (matchOrMismatch === COMPARE_MODES.COMPARE_MODE_MATCH ? true : false)
  );
};

export const getSinglePosColumn = ({
  pos,
  col,
  colorMode,
  refRow,
  compareMode,
  compareColor,
  colors,
  textColors,
}) => ({
  name: pos.toString(),
  key: col,
  cellClass: 'pos',
  sortable: true,
  width: 25,
  formatter: (val) => {
    // console.log(val);
    const row = val.row;
    let letter = row[col];
    let cellBgColor = 'transparent';
    let textColor = '#000';
    // Define the coloring behavior
    if (colorMode === COLOR_MODES.COLOR_MODE_COMPARE) {
      // Only color cells that meet the comparison condition (match/mismatch)
      // OR, if we're coloring by code, then always color the reference row
      if (
        (conditionCompare(row[col], refRow[col], compareMode) ||
          (row['group'] === GROUPS.REFERENCE_GROUP && colors !== null)) &&
        row[col] !== null
      ) {
        // If in dots mode, change letters, not colors
        if (compareColor === COMPARE_COLORS.COMPARE_COLOR_DOTS) {
          // Don't ever mask the reference with dots
          if (row['group'] !== GROUPS.REFERENCE_GROUP) {
            letter = '.';
          }
        } else if (colors === null) {
          cellBgColor = snapGeneHighlightColors[compareColor];
        } else {
          cellBgColor = colors[row[col]];
          textColor = textColors[row[col]];
        }
      }
    }
    // Just coloring by code
    else {
      cellBgColor = colors[row[col]];
    }

    return (
      <LetterCell value={letter} bgColor={cellBgColor} textColor={textColor} />
    );
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
  // {
  //   name: 'Show structure',
  //   key: null,
  //   sortable: false,
  //   width: 40,
  //   frozen: true,
  //   formatter: (val) => <AddToSidepanelCheckbox groupKey={val.row.group} />,
  // },
];
