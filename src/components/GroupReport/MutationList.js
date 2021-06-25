import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { format } from 'd3-format';

import { config } from '../../config';
import { getAllGenes, getAllProteins } from '../../utils/gene_protein';
import { reds, snpColorArray } from '../../constants/colors';
import { ASYNC_STATES, PLOT_DOWNLOAD_OPTIONS } from '../../constants/defs.json';

import EmptyPlot from '../Common/EmptyPlot';
import DropdownButton from '../Buttons/DropdownButton';
import SkeletonElement from '../Common/SkeletonElement';
import {
  MutationListContainer,
  MutationListHeader,
  OptionSelectContainer,
  OptionInputContainer,
  OptionCheckboxContainer,
  MutationListHeaderTable,
  MutationListHeaderEmpty,
  MutationListHeaderCell,
  MutationListTable,
  MutationRowContainer,
  MutationRowBar,
  MutationRowName,
  MutationRowHeatmapCellContainer,
} from './MutationList.styles';

const genes = getAllGenes();
const proteins = getAllProteins();
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

const MutationListRow = observer(
  ({
    segmentName,
    segmentColor,
    name,
    frequency,
    firstRow,
    numSnvsPerSegment,
    emptyRow,
  }) => {
    const { plotSettingsStore } = useStores();

    const toggleHiddenFeature = (featureName) => {
      plotSettingsStore.toggleReportMutationListHiddenItem(featureName);
    };

    const heatmapCells = [];
    if (!emptyRow) {
      // console.log(frequency.length);
      frequency.forEach((freq, i) => {
        heatmapCells.push(
          <MutationRowHeatmapCell
            key={`${name}-heatmap-cell-${i}`}
            freq={freq}
          />
        );
      });
    }

    return (
      <MutationRowContainer>
        {firstRow && (
          <MutationRowBar
            onClick={toggleHiddenFeature.bind(this, segmentName)}
            rowSpan={numSnvsPerSegment}
            barColor={segmentColor}
          >
            {segmentName}
          </MutationRowBar>
        )}
        <MutationRowName colSpan={emptyRow ? frequency.length + 1 : 1}>
          {name}
        </MutationRowName>
        {heatmapCells}
      </MutationRowContainer>
    );
  }
);
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
  const { groupDataStore, UIStore, plotSettingsStore } = useStores();

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

  // Use genes to group NT and gene_aa SNVs, and proteins for protein_aa SNVs
  const features =
    groupDataStore.groupSnvType === 'protein_aa' ? proteins : genes;

  features.forEach((feature, feature_i) => {
    // Get all SNVs for this gene, then sort by position/alt
    const groupFeatureSnvs = groupSnvFrequency
      .filter((groupSnv) => {
        if (groupDataStore.groupSnvType === 'dna') {
          // Include NT SNVs in this gene if it is contained in
          // any of the gene's NT segments
          // (Most genes will have one segment)
          return feature.segments.some(
            (featureNTRange) =>
              groupSnv.pos >= featureNTRange[0] &&
              groupSnv.pos <= featureNTRange[1]
          );
        } else if (groupDataStore.groupSnvType === 'gene_aa') {
          return groupSnv.gene === feature.name;
        } else if (groupDataStore.groupSnvType === 'protein_aa') {
          return groupSnv.protein === feature.name;
        }
      })
      .sort(sortByPosThenAlt);
    // console.log(feature.name, groupFeatureSnvs);

    // Make list of records to insert into master "matrix"
    const featureSnvRecords = groupFeatureSnvs
      .slice()
      // Unique SNVs
      .filter(
        (v, i, a) =>
          a.findIndex((element) => element.snv_name === v.snv_name) === i
      )
      .map((featureSnv) => {
        // Find fractional frequencies for each group
        const freqs = [];
        groupDataStore.selectedGroups.forEach((group) => {
          const matchingSnv = groupFeatureSnvs.find(
            (snv) => snv.snv_name === featureSnv.snv_name && snv.name === group
          );
          // 0 if the SNV record isn't found
          freqs.push(matchingSnv === undefined ? 0 : matchingSnv.fraction);
        });

        return {
          snv_name: featureSnv.snv_name,
          pos: featureSnv.pos,
          ref: featureSnv.ref,
          alt: featureSnv.alt,
          frequency: freqs,
        };
      })
      // Filter out SNVs that have all mutation frequencies below the threshold
      .filter((row) => {
        return row.frequency.some(
          (freq) => freq > plotSettingsStore.reportConsensusThreshold
        );
      });
    // console.log(feature.name, featureSnvRecords);

    // Push empty row for segments without SNVs
    if (featureSnvRecords.length === 0) {
      // Push nothing if we're hiding empty features
      if (plotSettingsStore.reportMutationListHideEmpty) {
        return;
      }

      rowItems.push(
        <MutationListRow
          key={`group-snv-empty-${feature.name}`}
          segmentName={feature.name}
          segmentColor={snpColorArray[feature_i % snpColorArray.length]}
          firstRow={true}
          frequency={new Array(groupDataStore.selectedGroups.length)}
          emptyRow={true}
        />
      );
      return;
    }
    // Push empty row for hidden features
    else if (
      plotSettingsStore.reportMutationListHidden.indexOf(feature.name) > -1
    ) {
      rowItems.push(
        <MutationListRow
          key={`group-snv-empty-${feature.name}`}
          segmentName={feature.name}
          segmentColor={snpColorArray[feature_i % snpColorArray.length]}
          firstRow={true}
          frequency={new Array(groupDataStore.selectedGroups.length)}
          emptyRow={true}
          name={`${featureSnvRecords.length} SNVs hidden...`}
        />
      );
      return;
    }

    featureSnvRecords.forEach((snv, i) => {
      const snvName =
        groupDataStore.groupSnvType === 'dna'
          ? snv.snv_name
          : snv.snv_name.split(':')[1];
      rowItems.push(
        <MutationListRow
          key={`group-snv-${feature.name}-${snv.snv_name}`}
          segmentName={feature.name}
          segmentColor={snpColorArray[feature_i % snpColorArray.length]}
          name={snvName}
          frequency={snv.frequency}
          firstRow={i === 0}
          numSnvsPerSegment={featureSnvRecords.length}
        />
      );
    });
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
  const { groupDataStore, plotSettingsStore } = useStores();

  // const onChangeActiveGroupType = (event) => {
  //   groupDataStore.updateActiveGroupType(event.target.value);
  // };
  const onChangeGroupSnvType = (event) => {
    groupDataStore.updateGroupSnvType(event.target.value);
  };
  const onChangeConsensusThreshold = (event) => {
    plotSettingsStore.setReportConsensusThreshold(event.target.value);
  };
  const onChangeHideEmpty = (event) => {
    plotSettingsStore.setReportMutationListHideEmpty(event.target.checked);
  };
  const handleDownloadSelect = (option) => {
    if (option === PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_DATA) {
      groupDataStore.downloadGroupSnvFrequencyData({
        group: groupDataStore.activeGroupType,
        snvType: groupDataStore.groupSnvType,
        consensusThreshold: 0,
      });
    }
  };

  if (groupDataStore.selectedGroups.length === 0) {
    return (
      <EmptyPlot height={250}>
        <p>No {groupDataStore.getActiveGroupTypePrettyName()}s selected</p>
      </EmptyPlot>
    );
  }

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
        {/* <OptionSelectContainer>
          <label>
            <select
              value={groupDataStore.activeGroupType}
              onChange={onChangeActiveGroupType}
            >
              {activeGroupTypeItems}
            </select>
          </label>
        </OptionSelectContainer> */}
        <OptionSelectContainer>
          <label>
            SNV Type
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
              value={plotSettingsStore.reportConsensusThreshold}
              onChange={onChangeConsensusThreshold}
              min={0}
              max={1}
              step={0.1}
            />
          </label>
        </OptionInputContainer>
        <div className="spacer"></div>
        <DropdownButton
          text={'Download'}
          options={[PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_DATA]}
          onSelect={handleDownloadSelect}
        />
      </MutationListHeader>
      <MutationListHeader>
        <OptionCheckboxContainer>
          <label>
            <input
              type="checkbox"
              name="mutation-list-hide-empty"
              checked={plotSettingsStore.reportMutationListHideEmpty}
              onChange={onChangeHideEmpty}
            />
            Hide{' '}
            {groupDataStore.groupSnvType === 'protein_aa'
              ? 'Proteins'
              : 'Genes'}{' '}
            without SNVs
          </label>
        </OptionCheckboxContainer>
      </MutationListHeader>
      <MutationListContent></MutationListContent>
    </MutationListContainer>
  );
});

export default MutationList;
