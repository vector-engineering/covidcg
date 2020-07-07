import React from 'react';
import styled from 'styled-components';
import PropTypes from 'prop-types';
import { snapGeneHighlightColors } from '../../utils/colors';
import { capitalize } from '../../utils/string';
import { connect, useStores } from '../../stores/connect';

import DropdownButton from '../Buttons/DropdownButton';
import SkeletonElement from '../SkeletonElement';
import { asyncStates } from '../../stores/uiStore';
import { observer } from 'mobx-react';

const ColorModeSelectLabel = styled.label`
  margin-right: 1em;
  select {
    margin-left: 0.65em;
    max-width: 200px;
    border-radius: 3px;
  }
`;

const ColorModeSelect = ({ dnaOrAa, colorMode, handleChange }) => {
  const colorModeElements = [
    <option key="compare" value="compare">
      Comparison to Reference
    </option>,
    <option key="code" value="code">
      {dnaOrAa === 'dna' ? '4' : '20'}-Color Code
    </option>,
  ];

  if (dnaOrAa === 'aa') {
    colorModeElements.push(
      <option key="clustal" value="clustal">
        ClustalX (Properties)
      </option>,
      <option key="zappo" value="zappo">
        Zappo (Physiochemical Properties)
      </option>,
      <option key="zhao-london" value="zhao-london">
        Zhao and London (Transmembrane Tendency)
      </option>
    );
  }

  return (
    <ColorModeSelectLabel>
      Color by:
      <select value={colorMode} onChange={handleChange}>
        {colorModeElements}
      </select>
    </ColorModeSelectLabel>
  );
};
ColorModeSelect.propTypes = {
  dnaOrAa: PropTypes.string.isRequired,
  colorMode: PropTypes.string.isRequired,
  handleChange: PropTypes.func.isRequired,
};

const CompareModeSelectLabel = styled.label`
  select {
    margin-left: 0.65em;
    max-width: 150px;
    border-radius: 3px;
  }
`;

const CompareModeSelect = ({
  disabled,
  dnaOrAa,
  compareMode,
  compareColor,
  handleModeChange,
  handleColorChange,
}) => {
  const colorOptions = Object.keys(snapGeneHighlightColors);
  const colorOptionElements = [];
  // Add the color-code option
  colorOptionElements.push(
    <option key="code" value="code">
      {dnaOrAa === 'dna' ? '4' : '20'}-color code
    </option>
  );

  if (dnaOrAa === 'aa') {
    colorOptionElements.push(
      <option key="clustal" value="clustal">
        ClustalX (Properties)
      </option>,
      <option key="zappo" value="zappo">
        Zappo (Physiochemical Properties)
      </option>,
      <option key="zhao-london" value="zhao-london">
        Zhao and London (Transmembrane Tendency)
      </option>
    );
  }

  colorOptions.forEach((color) => {
    colorOptionElements.push(
      <option key={color} value={color}>
        {capitalize(color)}
      </option>
    );
  });

  // Push dots option
  if (dnaOrAa === 'aa') {
    colorOptionElements.push(
      <option key="dots" value="dots">
        Dots
      </option>
    );
  }

  return (
    <CompareModeSelectLabel>
      Color bases that:
      <select
        disabled={disabled}
        value={compareMode}
        onChange={handleModeChange}
      >
        <option value="match">Match</option>
        <option value="mismatch">Don&apos;t Match</option>
      </select>
      <select
        disabled={disabled}
        value={compareColor}
        onChange={handleColorChange}
      >
        {colorOptionElements}
      </select>
    </CompareModeSelectLabel>
  );
};
CompareModeSelect.propTypes = {
  disabled: PropTypes.bool,
  dnaOrAa: PropTypes.string.isRequired,
  compareMode: PropTypes.string.isRequired,
  compareColor: PropTypes.string.isRequired,
  handleModeChange: PropTypes.func.isRequired,
  handleColorChange: PropTypes.func.isRequired,
};
CompareModeSelect.defaultProps = {
  disabled: false,
};

const DataTableOptions = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-start;

  margin-bottom: 10px;

  padding-left: 10px;
  padding-right: 24px;

  select {
    padding: 1px 4px;
    border-radius: 3px;
  }
`;

const Spacer = styled.div`
  flex-grow: 1;
`;

const SkeletonContainer = styled.div`
  display: flex;
  padding-top: 28px;
  padding-bottom: 28px;
  padding-left: 12px;
  padding-right: 24px;
`;

const TableOptions = observer(
  ({
    handleColorModeChange,
    handleCompareModeChange,
    handleCompareColorChange,
    colorMode,
    compareMode,
    compareColor,
  }) => {
    const { covidStore, uiStore } = useStores();

    const handleDownloadSelect = (option) => {
      if (option === 'Acknowledgements') {
        covidStore.downloadAcknowledgements();
      } else if (option === 'Aggregate Data') {
        covidStore.downloadAggCaseData();
      } else if (option === 'Seqeunces + Metadata') {
        covidStore.downloadSequencesAndMetadata();
      }
    };

    // if (
    //   uiStore.caseDataState === asyncStates.STARTED ||
    //   uiStore.aggCaseDataState === asyncStates.STARTED
    // ) {
    //   return (
    //     <SkeletonContainer>
    //       <SkeletonElement delay={2} height={'30px'} />
    //       &nbsp;&nbsp;
    //       <SkeletonElement delay={3} height={'30px'} />
    //       &nbsp;&nbsp;
    //       <SkeletonElement delay={4} height={'30px'} />
    //     </SkeletonContainer>
    //   );
    // }

    return (
      <DataTableOptions>
        <ColorModeSelect
          dnaOrAa={covidStore.dnaOrAa}
          colorMode={colorMode}
          handleChange={handleColorModeChange}
        />
        <CompareModeSelect
          disabled={colorMode === 'code'}
          dnaOrAa={covidStore.dnaOrAa}
          compareMode={compareMode}
          compareColor={compareColor}
          handleModeChange={handleCompareModeChange}
          handleColorChange={handleCompareColorChange}
        />
        <Spacer />
        <DropdownButton
          text={'Download'}
          options={[
            'Seqeunces + Metadata',
            'Aggregate Data',
            'Acknowledgements',
          ]}
          onSelect={handleDownloadSelect}
        />
      </DataTableOptions>
    );
  }
);

TableOptions.propTypes = {
  handleColorModeChange: PropTypes.func.isRequired,
  handleCompareModeChange: PropTypes.func.isRequired,
  handleCompareColorChange: PropTypes.func.isRequired,
  colorMode: PropTypes.string.isRequired,
  compareMode: PropTypes.string.isRequired,
  compareColor: PropTypes.string.isRequired,
};

export default connect(TableOptions);
