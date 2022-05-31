import React from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import {
  DNA_OR_AA,
  COORDINATE_MODES,
  GROUP_MUTATION,
  ASYNC_STATES,
} from '../../constants/defs.json';

import { formatMutation } from '../../utils/mutationUtils';
import { ISOToInt } from '../../utils/date';
import { getReferences } from '../../utils/reference';

import SkeletonElement from '../Common/SkeletonElement';
import DownloadDataButton from './DownloadDataButton';

import { Container, StatusText, Line, Sequence } from './StatusBox.styles';

const serializeCoordinates = (coordinateRanges) => {
  return coordinateRanges.map((coordRange) => coordRange.join('..')).join(', ');
};

const StatusBox = observer(() => {
  const { configStore, dataStore, UIStore } = useStores();

  let genomeSelection = '';
  const residuesOrBases =
    configStore.dnaOrAa === DNA_OR_AA.DNA ? 'Bases' : 'Residues';
  if (configStore.coordinateMode === COORDINATE_MODES.COORD_GENE) {
    genomeSelection = (
      <>
        Gene: <b>{configStore.selectedGene.name}</b>. {residuesOrBases}:{' '}
        <b>{serializeCoordinates(configStore.residueCoordinates)}</b>
      </>
    );
  } else if (configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
    genomeSelection = (
      <>
        Protein: <b>{configStore.selectedProtein.name}</b>. {residuesOrBases}:{' '}
        <b>{serializeCoordinates(configStore.residueCoordinates)}</b>
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
            : serializeCoordinates(configStore.getCoordinateRanges())}
        </b>
      </>
    );
  } else if (configStore.coordinateMode === COORDINATE_MODES.COORD_CUSTOM) {
    genomeSelection = (
      <>
        Custom coordinates:{' '}
        <b>{serializeCoordinates(configStore.customCoordinates)}</b>
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
        <b>{serializeCoordinates(configStore.getCoordinateRanges())}</b>
      </>
    );
  }

  const selectedGroupFields = [];
  Object.keys(configStore.selectedGroupFields).forEach((groupKey) => {
    if (configStore.selectedGroupFields[groupKey].length === 0) {
      return;
    }

    selectedGroupFields.push(
      <Line key={`status-box-selected-group-fields-${groupKey}`}>
        Selected {groupKey}s:{' '}
        {configStore.selectedGroupFields[groupKey].join(', ')}
      </Line>
    );
  });

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
