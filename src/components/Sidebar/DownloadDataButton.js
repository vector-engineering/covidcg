import React, { useState } from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import { ASYNC_STATES } from '../../constants/defs.json';
import { config } from '../../config';

import DropdownButton from '../Buttons/DropdownButton';
import SkeletonElement from '../Common/SkeletonElement';
import DownloadConsensusMutationsModal from '../Modals/DownloadConsensusMutationsModal';
import DownloadMetadataModal from '../Modals/DownloadMetadataModal';
import DownloadGenomesModal from '../Modals/DownloadGenomesModal';
import DownloadVariantTableModal from '../Modals/DownloadVariantTableModal';

import { ButtonContainer, DownloadButton } from './DownloadDataButton.styles';

const DOWNLOAD_OPTIONS = {
  AGGREGATE_DATA: 'Aggregate Data',
  GROUP_COUNTS: 'Group Counts',
  CONSENSUS_MUTATIONS: 'Consensus Mutations',
  VARIANT_TABLE: 'Variant Table',
  SELECTED_SEQUENCE_METADATA: 'Sequence Metadata',
  SELECTED_GENOMES: 'Selected Genomes',
};

const DownloadDataButton = observer(({ disabled, direction, style }) => {
  const { dataStore, UIStore } = useStores();

  const [activeModals, setActiveModals] = useState({
    downloadConsensusMutations: false,
    downloadMetadata: false,
    downloadGenomes: false,
    downloadVariantTable: false,
  });

  const showModal = (modal) => {
    setActiveModals({ ...activeModals, [modal]: true });
  };
  const hideModal = (modal) => {
    // Don't hide the modal until the download has finished
    if (
      UIStore.downloadState !== ASYNC_STATES.SUCCEEDED &&
      UIStore.downloadState !== ASYNC_STATES.FAILED
    ) {
      return;
    }
    setActiveModals({ ...activeModals, [modal]: false });
  };

  const handleDownloadSelect = (option) => {
    if (option === DOWNLOAD_OPTIONS.AGGREGATE_DATA) {
      dataStore.downloadAggSequences();
    } else if (option === DOWNLOAD_OPTIONS.GROUP_COUNTS) {
      dataStore.downloadAggGroup();
    } else if (option === DOWNLOAD_OPTIONS.CONSENSUS_MUTATIONS) {
      showModal('downloadConsensusMutations');
    } else if (option === DOWNLOAD_OPTIONS.VARIANT_TABLE) {
      showModal('downloadVariantTable');
    } else if (option === DOWNLOAD_OPTIONS.SELECTED_SEQUENCE_METADATA) {
      showModal('downloadMetadata');
    } else if (option === DOWNLOAD_OPTIONS.SELECTED_GENOMES) {
      showModal('downloadGenomes');
    }
  };

  const downloadOptions = [
    DOWNLOAD_OPTIONS.AGGREGATE_DATA,
    DOWNLOAD_OPTIONS.GROUP_COUNTS,
    DOWNLOAD_OPTIONS.CONSENSUS_MUTATIONS,
    DOWNLOAD_OPTIONS.VARIANT_TABLE,
  ];
  if (config.allow_metadata_download) {
    downloadOptions.push(DOWNLOAD_OPTIONS.SELECTED_SEQUENCE_METADATA);
  }
  if (config.allow_genome_download) {
    downloadOptions.push(DOWNLOAD_OPTIONS.SELECTED_GENOMES);
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
    <ButtonContainer style={style}>
      <DownloadConsensusMutationsModal
        isOpen={activeModals.downloadConsensusMutations}
        onRequestClose={hideModal.bind(this, 'downloadConsensusMutations')}
      />
      <DownloadMetadataModal
        isOpen={activeModals.downloadMetadata}
        onRequestClose={hideModal.bind(this, 'downloadMetadata')}
      />
      <DownloadGenomesModal
        isOpen={activeModals.downloadGenomes}
        onRequestClose={hideModal.bind(this, 'downloadGenomes')}
      />
      <DownloadVariantTableModal
        isOpen={activeModals.downloadVariantTable}
        onRequestClose={hideModal.bind(this, 'downloadVariantTable')}
      />
      <DropdownButton
        button={DownloadButton}
        text={'Download'}
        options={downloadOptions}
        onSelect={handleDownloadSelect}
        direction={direction}
        disabled={disabled}
      />
    </ButtonContainer>
  );
});

DownloadDataButton.propTypes = {
  disabled: PropTypes.bool,
  direction: PropTypes.string,
  style: PropTypes.object,
};
DownloadDataButton.defaultProps = {
  disabled: false,
  direction: 'left',
  style: {},
};

export default DownloadDataButton;
