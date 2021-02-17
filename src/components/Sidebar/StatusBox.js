import React from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import {
  DNA_OR_AA,
  COORDINATE_MODES,
  GROUP_SNV,
  ASYNC_STATES,
} from '../../constants/defs.json';
import { config } from '../../config';

import { formatSnv } from '../../utils/snpUtils';

import DropdownButton from '../Buttons/DropdownButton';
import SkeletonElement from '../Common/SkeletonElement';

import {
  Container,
  StatusText,
  ButtonContainer,
  Line,
  Sequence,
  DownloadButton,
} from './StatusBox.styles';

const serializeCoordinates = (coordinateRanges) => {
  return coordinateRanges.map((coordRange) => coordRange.join('..')).join(', ');
};

const DOWNLOAD_OPTIONS = {
  AGGREGATE_DATA: 'Aggregate Data',
  SELECTED_SEQUENCE_METADATA: 'Sequence Metadata',
};

const StatusBox = observer(() => {
  const { configStore, dataStore, UIStore } = useStores();

  const handleDownloadSelect = (option) => {
    if (option === DOWNLOAD_OPTIONS.AGGREGATE_DATA) {
      dataStore.downloadAggCaseData();
    } else if (option === DOWNLOAD_OPTIONS.SELECTED_SEQUENCE_METADATA) {
      dataStore.downloadSelectedSequenceMetadata();
    }
  };

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

  let selectedGroups = <b>None</b>;
  if (configStore.selectedGroups.length > 0) {
    if (configStore.groupKey === GROUP_SNV) {
      selectedGroups = configStore.selectedGroups
        .map((group) => {
          return (
            <b key={group.group}>
              {formatSnv(group.group, configStore.dnaOrAa)}
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

  const downloadOptions = [DOWNLOAD_OPTIONS.AGGREGATE_DATA];
  if (config.allow_metadata_download) {
    downloadOptions.push(DOWNLOAD_OPTIONS.SELECTED_SEQUENCE_METADATA);
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
      <StatusText>
        <Line>
          <b>{dataStore.numSequencesAfterAllFiltering}</b> sequences selected.
        </Line>
        <Line>
          Sequences grouped by <b>{configStore.getGroupLabel()}</b>.
        </Line>
        <Line>
          Viewing mutations on the{' '}
          <b>{configStore.dnaOrAa === DNA_OR_AA.DNA ? 'NT' : 'AA'}</b> level.
        </Line>
        <Line>
          Selected locations:{' '}
          <b>
            {configStore.selectedLocationNodes
              .map((node) => node.label)
              .join(', ')}
          </b>
        </Line>
        <Line>
          Date range: {configStore.startDate} – {configStore.endDate}
        </Line>
        <Line>Genome selection: {genomeSelection}</Line>
        <Line>
          Selected {configStore.getGroupLabel()}s: {selectedGroups}
        </Line>
      </StatusText>
      <ButtonContainer>
        <DropdownButton
          button={DownloadButton}
          text={'Download'}
          options={downloadOptions}
          onSelect={handleDownloadSelect}
          direction={'left'}
        />
      </ButtonContainer>
    </Container>
  );
});

export default StatusBox;
