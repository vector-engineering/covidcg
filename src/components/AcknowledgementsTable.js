import React from 'react';
import styled from 'styled-components';
import _ from 'underscore';

import { useStores } from '../stores/connect';

import DataGrid from 'react-data-grid';

import { getAckTextsFromAckIds } from '../utils/acknowledgements';

const AckContainer = styled.div``;

const AcknowledgementsTable = () => {
  const { covidStore } = useStores();

  let ackIds = _.pluck(covidStore.selectedRows, 'ack_id');
  // Get unique acknowledgement IDs
  ackIds = Array.from(new Set(ackIds));

  // Get the list of selected Accession IDs, and map to
  // acknowledgement texts
  let ackTexts = getAckTextsFromAckIds(ackIds);

  console.log(ackTexts);

  return (
    <AckContainer>
      <DataGrid
        columns={[
          {
            name: 'Originating Lab',
            key: 'Originating lab',
          },
          {
            name: 'Submitting Lab',
            key: 'Submitting lab',
          },
          {
            name: 'Authors',
            key: 'Authors',
          },
        ]}
        rows={ackTexts}
        height={200}
        headerRowHeight={45}
        filterRowHeight={45}
        rowHeight={25}
        minColumnWidth={25}
      />
    </AckContainer>
  );
};

export default AcknowledgementsTable;
