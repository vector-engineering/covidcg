import React, { useMemo, useState, useEffect } from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { format } from 'd3-format';

import { config } from '../../config';
import { getAllGenes, getAllProteins } from '../../utils/gene_protein';
import { reds, snpColorArray } from '../../constants/colors';
import { ASYNC_STATES } from '../../constants/defs.json';

import SkeletonElement from '../Common/SkeletonElement';
import {
  MutationListContainer,
  MutationListHeader,
  OptionSelectContainer,
  OptionInputContainer,
  MutationListHeaderTable,
  MutationListHeaderEmpty,
  MutationListHeaderCell,
  MutationListTable,
  MutationRowContainer,
  MutationRowBar,
  MutationRowName,
  MutationRowHeatmapCellContainer,
  MutationRowHeatmapEmptyCell,
} from './MutationList.styles';

const genes = getAllGenes();
const heatmapMin = 0.0;
const heatmapMax = 1.0;
const numColors = reds.length;
const countsFormatter = format('.0%');

const MutationRowHeatmapCell = ({ freq, percent }) => {
  // Find the color for this value
  let color = '#FFFFFF';
  // Add a bit extra since sometimes rounding errors can exclude the max value
  let interval = (heatmapMax - heatmapMin) / numColors + 0.00001;

  for (let i = 0; i < numColors; i++) {
    if (
      freq >= heatmapMin + i * interval &&
      freq <= heatmapMin + (i + 1) * interval
    ) {
      color = reds[i];
      break;
    }
  }

  // Don't show NaNs
  if (Number.isNaN(freq) || freq === null) {
    freq = '';
    color = 'transparent';
  }
  // Format percentages
  else if (percent === true) {
    freq = countsFormatter(freq);
  }

  return (
    <MutationRowHeatmapCellContainer bgColor={color}>
      {freq}
    </MutationRowHeatmapCellContainer>
  );
};
MutationRowHeatmapCell.propTypes = {
  freq: PropTypes.number.isRequired,
  percent: PropTypes.bool,
};
MutationRowHeatmapCell.defaultProps = {
  percent: true,
};

const MutationListRow = ({
  segmentName,
  segmentColor,
  name,
  frequency,
  firstRow,
  numSnvsPerSegment,
  emptyRow,
}) => {
  const heatmapCells = [];
  if (!emptyRow) {
    frequency.forEach((freq, i) => {
      heatmapCells.push(
        <MutationRowHeatmapCell key={`${name}-heatmap-cell-${i}`} freq={freq} />
      );
    });
  } else {
    console.log(frequency.length);
    heatmapCells.push(
      <MutationRowHeatmapEmptyCell
        key={`${segmentName}-empty-cell`}
        colSpan={frequency.length}
      />
    );
  }

  return (
    <MutationRowContainer>
      {firstRow && (
        <MutationRowBar rowSpan={numSnvsPerSegment} barColor={segmentColor}>
          {segmentName}
        </MutationRowBar>
      )}
      <MutationRowName>{name}</MutationRowName>
      {heatmapCells}
    </MutationRowContainer>
  );
};
MutationListRow.propTypes = {
  segmentName: PropTypes.string.isRequired,
  segmentColor: PropTypes.string.isRequired,
  name: PropTypes.string,
  frequency: PropTypes.arrayOf(PropTypes.number),
  firstRow: PropTypes.bool,
  numSnvsPerSegment: PropTypes.number,
  emptyRow: PropTypes.bool,
};
MutationListRow.defaultProps = {
  name: '',
  frequency: [],
  firstRow: false,
  numSnvsPerSegment: 1,
  emptyRow: false,
};

const MutationListContent = observer(() => {
  const { groupDataStore, UIStore } = useStores();

  // console.log(UIStore.groupSnvFrequencyState);

  if (UIStore.groupSnvFrequencyState !== ASYNC_STATES.SUCCEEDED) {
    return (
      <div
        style={{
          paddingTop: '12px',
          paddingRight: '24px',
          paddingLeft: '12px',
          paddingBottom: '24px',
        }}
      >
        <SkeletonElement delay={2} height={400} />
      </div>
    );
  }

  // Select group SNVs from the selected groups
  const groupSnvFrequency = groupDataStore.groupSnvFrequency[
    groupDataStore.activeGroupType
  ][groupDataStore.groupSnvType].filter((groupSnv) =>
    groupDataStore.selectedGroups.includes(groupSnv.name)
  );
  // console.log(groupSnvFrequency);

  // Restructure so that we have it in a matrix-ish format
  /*
  {
    // gene or protein
    S: [
      { 
        name: D614G
        ref, 
        alt,
        pos: 614
        frequency: [0.1, 0.9, 0.5] // fractional frequencies per group
      },
      ...
    ],
    ...
  ]
  */
  const rowItems = [];
  const sortByPosThenAlt = function (a, b) {
    if (a.pos === b.pos) {
      return a.alt > b.alt;
    } else {
      return a.pos - b.pos;
    }
  };
  genes.forEach((gene, gene_i) => {
    // Get all SNVs for this gene, then sort by position/alt
    const groupGeneSnvs = groupSnvFrequency
      .filter((groupSnv) => groupSnv.gene === gene.name)
      .sort(sortByPosThenAlt);
    // console.log(gene.name, groupGeneSnvs);

    // Make list of records to insert into master "matrix"
    const geneSnvRecords = groupGeneSnvs
      .slice()
      // Unique SNVs
      .filter(
        (v, i, a) =>
          a.findIndex((element) => element.snv_name === v.snv_name) === i
      )
      .map((geneSnv) => {
        // Find fractional frequencies for each group
        const freqs = [];
        groupDataStore.selectedGroups.forEach((group) => {
          const matchingSnv = groupGeneSnvs.find(
            (snv) => snv.snv_name === geneSnv.snv_name && snv.name === group
          );
          // 0 if the SNV record isn't found
          freqs.push(matchingSnv === undefined ? 0 : matchingSnv.fraction);
        });

        return {
          snv_name: geneSnv.snv_name,
          pos: geneSnv.pos,
          ref: geneSnv.ref,
          alt: geneSnv.alt,
          frequency: freqs,
        };
      });
    // console.log(gene.name, geneSnvRecords);

    geneSnvRecords.forEach((snv, i) => {
      rowItems.push(
        <MutationListRow
          key={`group-snv-${snv.snv_name}`}
          segmentName={gene.name}
          segmentColor={snpColorArray[gene_i % snpColorArray.length]}
          name={snv.snv_name.split(':')[1]}
          frequency={snv.frequency}
          firstRow={i === 0}
          numSnvsPerSegment={geneSnvRecords.length}
        />
      );
    });

    // Push empty row for segments without SNVs
    if (geneSnvRecords.length === 0) {
      rowItems.push(
        <MutationListRow
          key={`group-snv-empty-${gene.name}`}
          segmentName={gene.name}
          segmentColor={snpColorArray[gene_i % snpColorArray.length]}
          firstRow={true}
          frequency={new Array(groupDataStore.selectedGroups.length)}
          emptyRow={true}
        />
      );
    }
  });

  const headerItems = [];
  groupDataStore.selectedGroups.forEach((group) => {
    headerItems.push(
      <MutationListHeaderCell key={`mutation-list-table-head-${group}`}>
        <div>
          <span>{group}</span>
        </div>
      </MutationListHeaderCell>
    );
  });

  return (
    <>
      <MutationListHeaderTable ncols={groupDataStore.selectedGroups.length}>
        <thead>
          <tr>
            <MutationListHeaderEmpty colSpan={2} />
            {headerItems}
          </tr>
        </thead>
      </MutationListHeaderTable>
      <MutationListTable ncols={groupDataStore.selectedGroups.length}>
        <tbody>{rowItems}</tbody>
      </MutationListTable>
    </>
  );
});

const MutationList = observer(() => {
  const { groupDataStore } = useStores();

  const onChangeActiveGroupType = (event) => {
    groupDataStore.updateActiveGroupType(event.target.value);
  };
  const onChangeGroupSnvType = (event) => {
    groupDataStore.updateGroupSnvType(event.target.value);
  };
  const onChangeConsensusThreshold = (event) => {
    groupDataStore.updateConsensusThreshold(event.target.value);
  };

  const activeGroupTypeItems = [];
  Object.keys(config.group_cols).forEach((groupType) => {
    activeGroupTypeItems.push(
      <option key={`active-group-type-option-${groupType}`} value={groupType}>
        {config.group_cols[groupType].title}
      </option>
    );
  });

  return (
    <MutationListContainer>
      <MutationListHeader>
        <OptionSelectContainer>
          <label>
            <select
              value={groupDataStore.activeGroupType}
              onChange={onChangeActiveGroupType}
            >
              {activeGroupTypeItems}
            </select>
          </label>
        </OptionSelectContainer>
        <OptionSelectContainer>
          <label>
            <select
              value={groupDataStore.groupSnvType}
              onChange={onChangeGroupSnvType}
            >
              <option value={'dna'}>NT</option>
              <option value={'gene_aa'}>Gene AA</option>
              <option value={'protein_aa'}>Protein AA</option>
            </select>
          </label>
        </OptionSelectContainer>
        <OptionInputContainer>
          <label>
            Consensus Threshold
            <input
              type="number"
              value={groupDataStore.consensusThreshold}
              onChange={onChangeConsensusThreshold}
              min={0}
              max={1}
              step={0.1}
            />
          </label>
        </OptionInputContainer>
      </MutationListHeader>
      <MutationListContent></MutationListContent>
    </MutationListContainer>
  );
});

export default MutationList;
