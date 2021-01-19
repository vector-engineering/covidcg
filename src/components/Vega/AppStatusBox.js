import React from 'react';
import styled from 'styled-components';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import { ASYNC_STATES } from '../../constants/UI';
import { appConfig, DNA_OR_AA, COORDINATE_MODES, GROUP_SNV } from '../../constants/config';

import { formatSnv } from '../../utils/snpUtils';
import { intToISO } from '../../utils/date';

import DropdownButton from '../Buttons/DropdownButton';
import SkeletonElement from '../Common/SkeletonElement';

const Container = styled.div`
  display: grid;
  grid-template-columns: [col1] auto [col2] 110px [col3];
  grid-template-rows: [row1] auto [row2];

  margin: 0px 10px;
  padding: 5px 10px;
  border: 1px solid #CCC;
  border-radius: 5px;
`;

const LineColumn = styled.div`
  grid-row: row1 / row2;
  grid-column: col1 / col2;
`;

const ButtonColumn = styled.div`
  grid-row: row1 / row2;
  grid-column: col2 / col3;

  display: flex;
  flex-direction: column;
  align-items: flex-end;
`;

const Line = styled.p`
  font-size: 1em;
  font-weight: normal;
  margin: 0px;
`;

const Sequence = styled.span`
  font-family: monospace;
  display: inline;
  margin: 0px;
`;

const serializeCoordinates = (coordinateRanges) => {
  return coordinateRanges.map(coordRange => coordRange.join('..')).join(', ');
};

const DOWNLOAD_OPTIONS = {
  AGGREGATE_DATA: 'Aggregate Data',
  SELECTED_SEQUENCE_METADATA: 'Sequence Metadata'
};

const AppStatusBox = observer(() => {
  const { configStore, dataStore, UIStore } = useStores();

  const handleDownloadSelect = (option) => {
    if (option === DOWNLOAD_OPTIONS.AGGREGATE_DATA) {
      dataStore.downloadAggCaseData();
    } else if (option === DOWNLOAD_OPTIONS.SELECTED_SEQUENCE_METADATA) {
      dataStore.downloadSelectedSequenceMetadata();
    }
  };

  let genomeSelection = '';
  const residuesOrBases = configStore.dnaOrAa === DNA_OR_AA.DNA ? 'Bases' : 'Residues';
  if (configStore.coordinateMode === COORDINATE_MODES.COORD_GENE) {
    genomeSelection = (
    <>
      Gene: <b>{configStore.selectedGene.gene}</b>.{' '}
      {residuesOrBases}: <b>{serializeCoordinates(configStore.residueCoordinates)}</b>
    </>
    );
  } else if (configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
    genomeSelection = (
      <>
        Protein: <b>{configStore.selectedProtein.protein}</b>.{' '}
        {residuesOrBases}: <b>{serializeCoordinates(configStore.residueCoordinates)}</b>
      </>
    );
  } else if (configStore.coordinateMode === COORDINATE_MODES.COORD_PRIMER) {
    genomeSelection = (
      <>
        Primers:{' '}
        <b>
          {(
            configStore.selectedPrimers.length === 0 
              ? 'None' 
              : configStore.selectedPrimers.map(primer => primer.Name).join(', ')
          )}
        </b>.
        {' '}
        {configStore.selectedPrimers.length === 0 ? '' : (residuesOrBases + ': ')}
        <b>
          {(
            configStore.selectedPrimers.length === 0 
              ? '' 
              : serializeCoordinates(configStore.getCoordinateRanges())
          )}
        </b>
      </>
    );
  } else if (configStore.coordinateMode === COORDINATE_MODES.COORD_CUSTOM) {
    genomeSelection = (
      <>
        Custom coordinates: <b>{serializeCoordinates(configStore.customCoordinates)}</b>
      </>
    );
  } else if (configStore.coordinateMode === COORDINATE_MODES.COORD_SEQUENCE) {
    const sequenceList = configStore.customSequences.map((seq) => {
      return <Sequence key={seq}>{seq}</Sequence>
    }).reduce((prev, curr) => [prev, ', ', curr]);
    genomeSelection = (
      <>
        Matching sequence(s): {sequenceList}.{' '}
        {residuesOrBases}: <b>{serializeCoordinates(configStore.getCoordinateRanges())}</b>
      </>
    )
  }

  let dateRange = '';
  // Uninitialized date range is [-1, -1]
  if (configStore.dateRange[0] === -1 && configStore.dateRange[1] === -1) {
    dateRange = <b>All</b>;
  } else {
    dateRange = (
      <>
        <b>{intToISO(configStore.dateRange[0])}</b> – <b>{intToISO(configStore.dateRange[1])}</b>
      </>
    );
  }

  let selectedGroups = <b>None</b>;
  if (configStore.selectedGroups.length > 0) {
    if (configStore.groupKey === GROUP_SNV) {
      selectedGroups = configStore.selectedGroups.map((group) => {
        return <b key={group.group}>{formatSnv(group.group, configStore.dnaOrAa)}</b>;
      }).reduce((prev, curr) => [prev, ' + ', curr]);
    } else {
      selectedGroups = configStore.selectedGroups.map((group) => {
        return <b key={group.group}>{group.group}</b>;
      }).reduce((prev, curr) => [prev, ' + ', curr]);
    }
  }

  const downloadOptions = [DOWNLOAD_OPTIONS.AGGREGATE_DATA];
  if (appConfig.allow_metadata_download) {
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
        <SkeletonElement delay={2} height={100}>
        </SkeletonElement>
      </div>
    );
  }

  return (
    <Container>
      <LineColumn>
        <Line>
          <b>{dataStore.numSequencesAfterAllFiltering}</b> sequences selected. Sequences grouped by <b>{configStore.getGroupLabel()}</b>. Viewing mutations on the <b>{configStore.dnaOrAa === DNA_OR_AA.DNA ? 'NT' : 'AA'}</b> level.
        </Line>
        <Line>
          Selected locations: <b>{configStore.selectedLocationNodes.map((node) => node.label).join(', ')}</b>
        </Line>
        <Line>
          Date range: {dateRange}
        </Line>
        <Line>
          Genome selection: {genomeSelection}
        </Line>
        <Line>
          Selected {configStore.getGroupLabel()}s: {selectedGroups}
        </Line>
      </LineColumn>
      <ButtonColumn>
        <DropdownButton
          text={'Download'}
          options={downloadOptions}
          onSelect={handleDownloadSelect}
        />
      </ButtonColumn>
    </Container>
  );
});

export default AppStatusBox;
