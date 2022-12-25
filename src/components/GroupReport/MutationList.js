import React, { useState, useEffect } from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { format } from 'd3-format';

import { config } from '../../config';
import { reds, mutationColorArray } from '../../constants/colors';
import { ASYNC_STATES } from '../../constants/defs.json';
import { downloadBlobURL } from '../../utils/download';
import {
  buildFeatureMatrix,
  serializeFeatureMatrix,
} from './MutationList.utils';

import GroupSearch from './GroupSearch';
import EmptyPlot from '../Common/EmptyPlot';
import DropdownButton from '../Buttons/DropdownButton';
import SkeletonElement from '../Common/SkeletonElement';
import {
  MutationListContainer,
  MutationListHeader,
  MutationListTitle,
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
  MutationContentContainer,
  MutationInnerContainer,
  DeleteButton,
} from './MutationList.styles';
import { HelpButton, HelpText } from '../Common/Accordion.styles';

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
    numMutationsPerSegment,
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
            rowSpan={numMutationsPerSegment}
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
  numMutationsPerSegment: PropTypes.number,
  emptyRow: PropTypes.bool,
};
MutationListRow.defaultProps = {
  name: '',
  frequency: [],
  firstRow: false,
  numMutationsPerSegment: 1,
  emptyRow: false,
};

const DeleteButtonContainer = ({ group }) => {
  const { groupDataStore } = useStores();

  const onClick = () => {
    let selectedNodes = groupDataStore.selectedReportGroups;
    const index = selectedNodes.indexOf(group);
    if (index > -1) {
      selectedNodes.splice(index, 1);
    }

    // Update selectedNodes
    groupDataStore.applyPendingChanges({
      selectedReportGroups: selectedNodes.map((node) => {
        return node;
      }),
    });
  };

  return (
    <DeleteButton title="Deselect" onClick={onClick}>
      x
    </DeleteButton>
  );
};
DeleteButtonContainer.propTypes = {
  group: PropTypes.string.isRequired,
};

const MutationListContent = observer(({ featureMatrix }) => {
  const { groupDataStore, UIStore, plotSettingsStore } = useStores();

  if (UIStore.reportGroupMutationFrequencyState !== ASYNC_STATES.SUCCEEDED) {
    return (
      <div
        style={{
          paddingTop: '12px',
          paddingRight: '24px',
          paddingLeft: '12px',
          paddingBottom: '24px',
          overflow: 'auto',
        }}
      >
        <SkeletonElement delay={2} height={400} />
      </div>
    );
  }

  const rowItems = [];

  Object.keys(featureMatrix).forEach((featureName, feature_i) => {
    const featureMutationRecords = featureMatrix[featureName];

    // Push empty row for segments without mutations
    if (featureMutationRecords.length === 0) {
      // Push nothing if we're hiding empty features
      if (plotSettingsStore.reportMutationListHideEmpty) {
        return;
      }

      rowItems.push(
        <MutationListRow
          key={`group-mut-empty-${featureName}`}
          segmentName={featureName}
          segmentColor={
            mutationColorArray[feature_i % mutationColorArray.length]
          }
          firstRow={true}
          frequency={new Array(groupDataStore.selectedReportGroups.length)}
          emptyRow={true}
        />
      );
      return;
    } // Push empty row for hidden features
    else if (
      plotSettingsStore.reportMutationListHidden.indexOf(featureName) > -1
    ) {
      rowItems.push(
        <MutationListRow
          key={`group-mut-empty-${featureName}`}
          segmentName={featureName}
          segmentColor={
            mutationColorArray[feature_i % mutationColorArray.length]
          }
          firstRow={true}
          frequency={new Array(groupDataStore.selectedReportGroups.length)}
          emptyRow={true}
          name={`${featureMutationRecords.length} mutations hidden...`}
        />
      );
      return;
    }

    featureMutationRecords.forEach((mut, i) => {
      const mutName =
        groupDataStore.reportGroupMutationType === 'dna'
          ? mut.mutation_name
          : mut.mutation_name.split(':')[1];

      rowItems.push(
        <MutationListRow
          key={`group-mut-${featureName}-${mut.mutation_name}`}
          segmentName={featureName}
          segmentColor={
            mutationColorArray[feature_i % mutationColorArray.length]
          }
          name={mutName}
          frequency={mut.frequency}
          firstRow={i === 0}
          numMutationsPerSegment={featureMutationRecords.length}
        />
      );
    });
  });

  const headerItems = [];
  const deleteButtons = [];
  groupDataStore.selectedReportGroups.forEach((group, i) => {
    headerItems.push(
      <MutationListHeaderCell key={`mutation-list-table-head-${group}`}>
        <div>
          <span>{group}</span>
        </div>
      </MutationListHeaderCell>
    );
    deleteButtons.push(
      <DeleteButtonContainer key={`delete-${group}-${i}`} group={group} />
    );
  });

  return (
    <MutationContentContainer>
      <MutationListHeaderTable
        ncols={groupDataStore.selectedReportGroups.length}
      >
        <thead>
          <tr>
            <GroupSearch />
            {headerItems}
          </tr>
          <tr>
            <MutationListHeaderEmpty colSpan={2} />
            {deleteButtons}
          </tr>
        </thead>
      </MutationListHeaderTable>
      <MutationListTable ncols={groupDataStore.selectedReportGroups.length}>
        <tbody>{rowItems}</tbody>
      </MutationListTable>
    </MutationContentContainer>
  );
});
MutationListContent.propTypes = {
  featureMatrix: PropTypes.object.isRequired,
};

const DOWNLOAD_OPTIONS = {
  MUTATION_DATA: 'Mutation Data',
  MUTATION_DATA_MATRIX: 'Mutation Data (Matrix)',
};

const MutationList = observer(() => {
  const { groupDataStore, plotSettingsStore, UIStore } = useStores();

  const [state, setState] = useState({
    showHelp: false,
    featureMatrix: buildFeatureMatrix({
      reportActiveReference: groupDataStore.reportActiveReference,
      reportGroupMutationFrequency: groupDataStore.reportGroupMutationFrequency,
      activeReportGroupType: groupDataStore.activeReportGroupType,
      reportGroupMutationType: groupDataStore.reportGroupMutationType,
      selectedReportGroups: groupDataStore.selectedReportGroups,
      reportConsensusThreshold: plotSettingsStore.reportConsensusThreshold,
    }),
  });

  // const onChangeActiveReportGroupType = (event) => {
  //   groupDataStore.applyPendingChanges({ activeReportGroupType: event.target.value });
  // };
  const onChangeReportGroupMutationType = (event) => {
    groupDataStore.applyPendingChanges({
      reportGroupMutationType: event.target.value,
    });
  };
  const onChangeConsensusThreshold = (event) => {
    plotSettingsStore.setReportConsensusThreshold(event.target.value);
  };
  const onChangeHideEmpty = (event) => {
    plotSettingsStore.setReportMutationListHideEmpty(event.target.checked);
  };
  const handleDownloadSelect = (option) => {
    if (option === DOWNLOAD_OPTIONS.MUTATION_DATA) {
      groupDataStore.downloadGroupMutationFrequencyData({
        group: groupDataStore.activeReportGroupType,
        mutationType: groupDataStore.reportGroupMutationType,
        consensusThreshold: 0,
        selectedReference: groupDataStore.reportActiveReference,
      });
    } else if (option === DOWNLOAD_OPTIONS.MUTATION_DATA_MATRIX) {
      const blob = new Blob([
        serializeFeatureMatrix({
          featureMatrix: state.featureMatrix,
          selectedReportGroups: groupDataStore.selectedReportGroups,
        }),
      ]);
      const url = URL.createObjectURL(blob);
      downloadBlobURL(url, 'mutation_data_matrix.csv');
    }
  };
  const toggleHelp = (e) => {
    e.preventDefault();
    setState({ ...state, showHelp: !state.showHelp });
  };

  useEffect(() => {
    setState({
      ...state,
      featureMatrix: buildFeatureMatrix({
        reportActiveReference: groupDataStore.reportActiveReference,
        reportGroupMutationFrequency:
          groupDataStore.reportGroupMutationFrequency,
        activeReportGroupType: groupDataStore.activeReportGroupType,
        reportGroupMutationType: groupDataStore.reportGroupMutationType,
        selectedReportGroups: groupDataStore.selectedReportGroups,
        reportConsensusThreshold: plotSettingsStore.reportConsensusThreshold,
      }),
    });
  }, [
    UIStore.reportGroupMutationFrequencyState,
    groupDataStore.selectedReportGroups,
  ]);

  if (groupDataStore.selectedReportGroups.length === 0) {
    return (
      <EmptyPlot height={250}>
        <p>
          No {groupDataStore.getActiveReportGroupTypePrettyName()}s selected
        </p>
      </EmptyPlot>
    );
  }

  const activeReportGroupTypeItems = [];
  Object.keys(config.group_cols).forEach((groupType) => {
    activeReportGroupTypeItems.push(
      <option key={`active-group-type-option-${groupType}`} value={groupType}>
        {config.group_cols[groupType].title}
      </option>
    );
  });

  return (
    <MutationListContainer>
      <MutationListHeader>
        <MutationListTitle>Lineage Mutations</MutationListTitle>
        <HelpButton onClick={toggleHelp}>Show Help</HelpButton>
      </MutationListHeader>
      <MutationListHeader>
        <HelpText show={state.showHelp}>
          <ul>
            <li>
              Use &quot;Mutation Type&quot; to toggle between nucleotide and
              amino acid mutation formats
            </li>
            <li>
              &quot;Consensus Threshold&quot; hides low-prevalence mutations.
              Mutations with less than this fraction of prevalence in <i>all</i>{' '}
              selected {groupDataStore.getReportGroupMutationTypePrettyName()}s
              will be filtered out.
            </li>
            <li>
              Note: We define ORF1a and ORF1ab as separate genes. In
              &quot;NT&quot; or &quot;Gene AA&quot; mode, a mutation in ORF1a
              will also be listed in ORF1ab. By default, ORF1a is hidden to
              avoid this confusion.
            </li>
            <li>
              Switch to &quot;Protein AA&quot; mode to see mutations in the
              context of proteins (i.e., NSPs)
            </li>
          </ul>
        </HelpText>
      </MutationListHeader>
      <MutationListHeader>
        {/* <OptionSelectContainer>
          <label>
            <select
              value={groupDataStore.activeReportGroupType}
              onChange={onChangeActiveReportGroupType}
            >
              {activeReportGroupTypeItems}
            </select>
          </label>
        </OptionSelectContainer> */}
        <OptionSelectContainer>
          <label>
            Mutation Type
            <select
              value={groupDataStore.reportGroupMutationType}
              onChange={onChangeReportGroupMutationType}
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
          options={[
            DOWNLOAD_OPTIONS.MUTATION_DATA,
            DOWNLOAD_OPTIONS.MUTATION_DATA_MATRIX,
          ]}
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
            {groupDataStore.reportGroupMutationType === 'protein_aa'
              ? 'Proteins'
              : 'Genes'}{' '}
            without mutations
          </label>
        </OptionCheckboxContainer>
      </MutationListHeader>
      <MutationInnerContainer>
        <MutationListContent featureMatrix={state.featureMatrix} />
      </MutationInnerContainer>
    </MutationListContainer>
  );
});

export default MutationList;
