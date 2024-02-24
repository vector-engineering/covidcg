import React from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { config } from '../../config';

import {
  DNA_OR_AA,
  COORDINATE_MODES,
  GROUP_MUTATION,
  ASYNC_STATES,
} from '../../constants/defs.json';

import { coordsToText } from '../../utils/coordinates';
import { formatMutation } from '../../utils/mutationUtils';
import { ISOToInt } from '../../utils/date';
import { getReferences } from '../../utils/reference';

import SkeletonElement from '../Common/SkeletonElement';
import DownloadDataButton from './DownloadDataButton';

import { Container, StatusText, Line, Sequence } from './StatusBox.styles';

const serializeNTCoordinates = (coordinateRanges) => {
  return coordinateRanges
    .map((coordRange) => `${coordRange[0]}:${coordRange[1]}..${coordRange[2]}`)
    .join(', ');
};
const serializeAACoordinates = (coordinateRanges) => {
  return coordinateRanges.map((coordRange) => coordRange.join('..')).join(', ');
};

const StatusBox = observer(() => {
  const { configStore, dataStore, UIStore } = useStores();

  let qualityFilters = '';
  const qualityFields = ['length', 'percent_ambiguous'];
  const qualityFieldKeys = {
    length: 'sequenceLengthRange',
    percent_ambiguous: 'percentAmbiguousRange',
  };
  const qualityFieldNames = {
    length: 'Sequence Length',
    percent_ambiguous: '% Ambiguous (N) Bases',
  };
  // Only render this if we have the quality filters available
  if (config['virus'] === 'sars2') {
    const qualityFilterItems = [];
    qualityFields.forEach((field) => {
      const rng = configStore[qualityFieldKeys[field]];
      let suffix = ' bases';
      if (field === 'percent_ambiguous') {
        suffix = '%';
      }
      // If first value of the range undefined, then only maximum is set
      if (rng[0] === null) {
        qualityFilterItems.push(
          <Line key={`quality-status-${field}`}>
            {qualityFieldNames[field]} ≤{' '}
            <b>
              {rng[1]}
              {suffix}
            </b>
          </Line>
        );
      }
      // If second value of the range undefined, then only minimum is set
      else if (rng[1] === null) {
        qualityFilterItems.push(
          <Line key={`quality-status-${field}`}>
            {qualityFieldNames[field]} ≥{' '}
            <b>
              {rng[0]}
              {suffix}
            </b>
          </Line>
        );
      }
      // If both values of the range are defined, then both minimum and maximum are set
      else {
        qualityFilterItems.push(
          <Line key={`quality-status-${field}`}>
            {qualityFieldNames[field]}:{' '}
            <b>
              {rng[0]}
              {suffix}
            </b>{' '}
            –{' '}
            <b>
              {rng[1]}
              {suffix}
            </b>
          </Line>
        );
      }
    });
    qualityFilters = <>{qualityFilterItems}</>;
  }

  const selectedGroupFields = [];
  Object.keys(configStore.selectedGroupFields).forEach((groupKey) => {
    if (configStore.selectedGroupFields[groupKey].length === 0) {
      return;
    }

    const selectedGroupFieldItems = [];
    configStore.selectedGroupFields[groupKey].forEach((group, i) => {
      selectedGroupFieldItems.push(
        <b key={`status-selected-group-${groupKey}-${group}}`}>{group}</b>
      );
      if (i < configStore.selectedGroupFields[groupKey].length - 1) {
        selectedGroupFieldItems.push(',');
      }
    });

    selectedGroupFields.push(
      <Line key={`status-box-selected-group-fields-${groupKey}`}>
        Selected {groupKey}s: {selectedGroupFieldItems}
      </Line>
    );
  });

  let genomeSelection = '';
  const residuesOrBases =
    configStore.dnaOrAa === DNA_OR_AA.DNA ? 'Bases' : 'Residues';
  const coordRange =
    configStore.dnaOrAa === DNA_OR_AA.DNA
      ? serializeNTCoordinates(configStore.getCoordinateRanges())
      : serializeAACoordinates(configStore.residueCoordinates);
  if (configStore.coordinateMode === COORDINATE_MODES.COORD_GENE) {
    genomeSelection = (
      <>
        Gene: <b>{configStore.selectedGene.name}</b>. {residuesOrBases}:{' '}
        <b>{coordRange}</b>
      </>
    );
  } else if (configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
    genomeSelection = (
      <>
        Protein: <b>{configStore.selectedProtein.name}</b>. {residuesOrBases}:{' '}
        <b>{coordRange}</b>
      </>
    );
  } else if (configStore.coordinateMode === COORDINATE_MODES.COORD_PRIMER) {
    genomeSelection = (
      <>
        Primers:{' '}
        <b>
          {configStore.selectedPrimers.length === 0
            ? 'None'
            : configStore.selectedPrimers
                .map((primer) => primer.Name)
                .join(', ')}
        </b>
        .{' '}
        {configStore.selectedPrimers.length === 0 ? '' : residuesOrBases + ': '}
        <b>
          {configStore.selectedPrimers.length === 0
            ? ''
            : serializeNTCoordinates(configStore.getCoordinateRanges())}
        </b>
      </>
    );
  } else if (configStore.coordinateMode === COORDINATE_MODES.COORD_CUSTOM) {
    genomeSelection = (
      <>
        Custom coordinates: <b>{coordsToText(configStore.customCoordinates)}</b>
      </>
    );
  } else if (configStore.coordinateMode === COORDINATE_MODES.COORD_SEQUENCE) {
    const sequenceList = configStore.customSequences
      .map((seq) => {
        return <Sequence key={seq}>{seq}</Sequence>;
      })
      .reduce((prev, curr) => [prev, ', ', curr]);
    genomeSelection = (
      <>
        Matching sequence(s): {sequenceList}. {residuesOrBases}:{' '}
        <b>{serializeNTCoordinates(configStore.getCoordinateRanges())}</b>
      </>
    );
  }

  let selectedGroups = <b>None</b>;
  if (configStore.selectedGroups.length > 0) {
    if (configStore.groupKey === GROUP_MUTATION) {
      selectedGroups = configStore.selectedGroups
        .map((group) => {
          return (
            <b key={group.group}>
              {formatMutation(group.group, configStore.dnaOrAa)}
            </b>
          );
        })
        .reduce((prev, curr) => [prev, ' + ', curr]);
    } else {
      selectedGroups = configStore.selectedGroups
        .map((group) => {
          return <b key={group.group}>{group.group}</b>;
        })
        .reduce((prev, curr) => [prev, ' + ', curr]);
    }
  }

  if (UIStore.caseDataState === ASYNC_STATES.STARTED) {
    return (
      <div
        style={{
          paddingTop: '12px',
          paddingRight: '24px',
          paddingLeft: '12px',
          paddingBottom: '24px',
        }}
      >
        <SkeletonElement delay={2} height={100}></SkeletonElement>
      </div>
    );
  }

  return (
    <Container>
      <DownloadDataButton />
      <StatusText>
        <Line>
          <b>{dataStore.numSequencesAfterAllFiltering}</b> sequences fetched in{' '}
          {dataStore.timeToFetch} s.
        </Line>
        <Line>
          Sequences grouped by <b>{configStore.getGroupLabel()}</b>.
        </Line>
        <Line>
          Reference genome: <b>{configStore.selectedReference}</b> (
          {getReferences()[configStore.selectedReference]['description']}).
        </Line>
        {qualityFilters}
        <Line>
          <b>{configStore.selectedLocationNodes.length}</b> selected locations:{' '}
          <b>
            {configStore.selectedLocationNodes
              .map((node) => node.label)
              .join(', ')}
          </b>
        </Line>
        <Line>
          Date range: <b>{configStore.startDate}</b> – 
          <b>{configStore.endDate}</b> (
          {(ISOToInt(configStore.endDate) - ISOToInt(configStore.startDate)) /
            (1000 * 60 * 60 * 24)}{' '}
          days)
        </Line>
        {selectedGroupFields}
        {configStore.groupKey === GROUP_MUTATION && (
          <Line>Genome selection: {genomeSelection}</Line>
        )}
        <Line>
          Selected {configStore.getGroupLabel()}s: {selectedGroups}
        </Line>
      </StatusText>
    </Container>
  );
});

export default StatusBox;
