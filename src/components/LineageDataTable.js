/* eslint-disable react/display-name */
import React, { useState } from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import _ from 'underscore';

import DataTable from 'react-data-table-component';
import HeatmapCell from './Cells/HeatmapCell';
import LetterCell from './Cells/LetterCell';
import Button from './Button';
import AddToSidepanelCheckbox from './AddToSidepanelCheckbox';
import { useStores } from '../stores/connect';
import { nanmin, nanmax } from '../utils/math';
import {
  snapGeneHighlightColors,
  snapGeneNTColors,
  shingAAColors,
} from '../utils/colors';
import { capitalize } from '../utils/string';

const DataTableContainer = styled.div`
  span.position-title {
    font-size: 12px;
    font-weight: 500;
  }
`;

const DataTableOptions = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-start;

  margin-bottom: 10px;

  padding-left: 10px;
  padding-right: 10px;

  select {
    padding: 1px 4px;
  }
`;

const ColorModeSelectLabel = styled.label`
  margin-right: 1em;
  select {
    margin-left: 0.65em;
  }
`;

const CompareModeSelectLabel = styled.label`
  select {
    margin-left: 0.65em;
  }
`;

const Spacer = styled.div`
  flex-grow: 1;
`;

const StyledDataTable = styled(DataTable)`
  margin-top: 3px;

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

  .rdt_TableCol_Sortable[id^='column-pos_'] {
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

      & > div {
        overflow: visible;
      }
    }
  }
`;

const ColorModeSelect = ({ dnaOrAa, colorMode, handleChange }) => {
  const colorModeElements = [
    <option key="compare" value="compare">
      Comparison to Reference
    </option>,
    <option key="code" value="code">
      {dnaOrAa === 'dna' ? '4' : '20'}-Color Code
    </option>,
  ];

  return (
    <ColorModeSelectLabel>
      Color by:
      <select value={colorMode} onChange={handleChange}>
        {colorModeElements}
      </select>
    </ColorModeSelectLabel>
  );
};
ColorModeSelect.propTypes = {
  dnaOrAa: PropTypes.string.isRequired,
  colorMode: PropTypes.string.isRequired,
  handleChange: PropTypes.func.isRequired,
};

const CompareModeSelect = ({
  disabled,
  dnaOrAa,
  compareMode,
  compareColor,
  handleModeChange,
  handleColorChange,
}) => {
  const colorOptions = Object.keys(snapGeneHighlightColors);
  const colorOptionElements = [];
  // Add the color-code option
  colorOptionElements.push(
    <option key="code" value="code">
      {dnaOrAa === 'dna' ? '4' : '20'}-color code
    </option>
  );
  colorOptions.forEach((color) => {
    colorOptionElements.push(
      <option key={color} value={color}>
        {capitalize(color)}
      </option>
    );
  });

  return (
    <CompareModeSelectLabel>
      Color bases that:
      <select
        disabled={disabled}
        value={compareMode}
        onChange={handleModeChange}
      >
        <option value="match">Match</option>
        <option value="mismatch">Don&apos;t Match</option>
      </select>
      <select
        disabled={disabled}
        value={compareColor}
        onChange={handleColorChange}
      >
        {colorOptionElements}
      </select>
    </CompareModeSelectLabel>
  );
};
CompareModeSelect.propTypes = {
  disabled: PropTypes.bool,
  dnaOrAa: PropTypes.string.isRequired,
  compareMode: PropTypes.string.isRequired,
  compareColor: PropTypes.string.isRequired,
  handleModeChange: PropTypes.func.isRequired,
  handleColorChange: PropTypes.func.isRequired,
};
CompareModeSelect.defaultProps = {
  disabled: false,
};

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

const conditionCompare = (base, refBase, matchOrMismatch) => {
  // Flip the XOR (XNOR)
  return !((base === refBase) ^ (matchOrMismatch === 'match' ? true : false));
};

const LineageDataTable = observer(() => {
  const { covidStore } = useStores();
  // Color by 'compare': Comparison to reference, or 'code': With a defined color code
  const [colorMode, setColorMode] = useState('compare');
  // 'match' or 'mismatch'
  const [compareMode, setCompareMode] = useState('mismatch');
  const [compareColor, setCompareColor] = useState('yellow');

  const handleColorModeChange = (event) => setColorMode(event.target.value);
  const handleCompareModeChange = (event) => setCompareMode(event.target.value);
  const handleCompareColorChange = (event) =>
    setCompareColor(event.target.value);

  const downloadAcknowledgements = () => covidStore.downloadAcknowledgements();

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

  // The offset of the position title will depend on what columns
  // are present
  let positionTitleOffset = 0;

  // Build initial columns
  const columns = [];

  // For lineage grouping, add lineage column
  if (covidStore.groupKey === 'lineage') {
    columns.push({
      name: 'Lineage',
      selector: 'group',
      sortable: true,
      width: '100px',
      style: { fontWeight: '700' },
    });
    positionTitleOffset += 100;
  }

  // For SNP grouping, add each SNP chunk as its own column
  if (covidStore.groupKey === 'snp') {
    // Add the gene column, if we're in AA mode
    if (covidStore.dnaOrAa === 'aa') {
      columns.push({
        name: 'Gene',
        selector: 'gene',
        sortable: true,
        width: '50px',
        style: { fontWeight: '700' },
      });
      positionTitleOffset += 50;
    }
    // Add the position column
    // We don't need as much space for this, for AA mode
    columns.push({
      name: covidStore.dnaOrAa === 'dna' ? 'Position' : 'Index',
      selector: 'pos',
      sortable: true,
      width: covidStore.dnaOrAa === 'dna' ? '60px' : '45px',
      style: { fontWeight: '700' },
    });
    positionTitleOffset += covidStore.dnaOrAa === 'dna' ? 60 : 45;
    // Add the Ref and Alt columns
    columns.push({
      name: 'Ref',
      selector: 'ref',
      sortable: true,
      width: '24px',
      cell: (row) => <LetterCell value={row['ref']} />,
    });
    columns.push({
      name: 'Alt',
      selector: 'alt',
      sortable: true,
      width: '24px',
      cell: (row) => <LetterCell value={row['alt']} />,
    });
    positionTitleOffset += 48;
  }

  // Push the quantitative columns
  columns.push({
    name: 'Seqs',
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
  });
  positionTitleOffset += 60;
  columns.push({
    name: '% Seqs',
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
  });
  positionTitleOffset += 75;

  // Push sidepanel column
  columns.push({
    name: 'Show structure',
    selector: null,
    sortable: false,
    width: '80px',
    cell: (row) => {
      return <AddToSidepanelCheckbox groupKey={row.group} />;
    },
  });
  positionTitleOffset += 80;

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

    columns.push({
      name: pos.toString(),
      selector: col,
      sortable: false,
      width: '24px',
      cell: (row) => {
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
      center: true,
      compact: true,
    });
  });

  // Default sorting behavior
  let defaultSortBy;
  // For lineage, sort by lineage
  if (covidStore.groupKey === 'lineage') {
    defaultSortBy = 'group';
  }
  // For SNPs, sort by position
  else if (covidStore.groupKey === 'snp') {
    defaultSortBy = 'pos';
  }

  return (
    <DataTableContainer>
      <DataTableOptions>
        <ColorModeSelect
          dnaOrAa={covidStore.dnaOrAa}
          colorMode={colorMode}
          handleChange={handleColorModeChange}
        />
        <CompareModeSelect
          disabled={colorMode === 'code'}
          dnaOrAa={covidStore.dnaOrAa}
          compareMode={compareMode}
          compareColor={compareColor}
          handleModeChange={handleCompareModeChange}
          handleColorChange={handleCompareColorChange}
        />
        <Spacer />
        <Button onClick={downloadAcknowledgements}>
          Download Acknowledgements
        </Button>
      </DataTableOptions>
      <span
        className="position-title"
        style={{ marginLeft: positionTitleOffset }}
      >
        {covidStore.dnaOrAa === 'dna' ? 'Genomic Coordinate' : 'Residue Index'}
      </span>
      <StyledDataTable
        className="data-table"
        data={covidStore.caseDataAggGroup}
        keyField="group"
        columns={columns}
        striped={true}
        highlightOnHover={true}
        dense={true}
        // fixedHeader={true}
        // fixedHeaderScrollHeight={'400px'}
        pagination={false}
        defaultSortField={defaultSortBy}
        defaultSortAsc={true}
        conditionalRowStyles={[
          {
            when: (row) => row.group == 'Reference',
            style: { backgroundColor: '#dff3fe !important;' },
          },
        ]}
        sortFunction={sortFn}
      />
    </DataTableContainer>
  );
});

LineageDataTable.displayName = 'LineageDataTable';

export default LineageDataTable;
