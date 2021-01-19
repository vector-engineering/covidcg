import React from 'react';
import styled from 'styled-components';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import { DNA_OR_AA, COORDINATE_MODES, GROUP_SNV } from '../../constants/config';

import { formatSnv } from '../../utils/snpUtils';

const Container = styled.div`
  margin: 0px 10px;
  padding: 5px 10px;
  border: 1px solid #CCC;
  border-radius: 5px;
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

// I can't believe I have to write my own date formatting code. 
// wtf happened to standard C formatting
const dateToISO = (dateInUnixMS) => {
  const dateObj = new Date(dateInUnixMS);
  return (
    dateObj.getFullYear() + '-' +
    // Add 1 since month is in [0, 11], then pad to 2 characters
    (dateObj.getMonth() + 1).toString().padStart(2, '0') + '-' +
    // Pad to 2 characters. Date is [1, 31].
    dateObj.getDate().toString().padStart(2, '0')
  );
}

const AppStatusBox = observer(() => {
  const { configStore, dataStore } = useStores();

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
        <b>{dateToISO(configStore.dateRange[0])}</b> – <b>{dateToISO(configStore.dateRange[1])}</b>
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

  return (
    <Container>
      <Line>
        <b>{dataStore.filteredCaseData.length}</b> sequences selected. Sequences grouped by <b>{configStore.getGroupLabel()}</b>. Viewing mutations on the <b>{configStore.dnaOrAa === DNA_OR_AA.DNA ? 'NT' : 'AA'}</b> level.
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
    </Container>
  );
});

export default AppStatusBox;
