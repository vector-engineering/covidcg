import React from 'react';
import styled from 'styled-components';
import { snapGeneHighlightColors } from '../../constants/colors';
import { capitalize } from '../../utils/string';
import { useStores } from '../../stores/connect';
import {
  COLOR_MODES,
  COMPARE_MODES,
  COMPARE_COLORS,
} from '../../stores/plotSettingsStore';

import DropdownButton from '../Buttons/DropdownButton';
import { observer } from 'mobx-react';

const ColorModeSelectLabel = styled.label`
  margin-right: 1em;
  select {
    margin-left: 0.65em;
    max-width: 200px;
    border-radius: 3px;
  }
`;

const ColorModeSelect = observer(() => {
  const { plotSettingsStore, configStore } = useStores();

  const handleColorModeChange = (event) =>
    plotSettingsStore.setTableColorMode(event.target.value);

  const colorModeElements = [
    <option
      key={COLOR_MODES.COLOR_MODE_COMPARE}
      value={COLOR_MODES.COLOR_MODE_COMPARE}
    >
      Comparison to Reference
    </option>,
    <option
      key={COLOR_MODES.COLOR_MODE_CODE}
      value={COLOR_MODES.COLOR_MODE_CODE}
    >
      {configStore.dnaOrAa === 'dna' ? '4' : '20'}-Color Code
    </option>,
  ];

  if (configStore.dnaOrAa === 'aa') {
    colorModeElements.push(
      <option
        key={COLOR_MODES.COLOR_MODE_CLUSTAL}
        value={COLOR_MODES.COLOR_MODE_CLUSTAL}
      >
        ClustalX (Properties)
      </option>,
      <option
        key={COLOR_MODES.COLOR_MODE_ZAPPO}
        value={COLOR_MODES.COLOR_MODE_ZAPPO}
      >
        Zappo (Physiochemical Properties)
      </option>,
      <option
        key={COLOR_MODES.COLOR_MODE_ZHAO_LONDON}
        value={COLOR_MODES.COLOR_MODE_ZHAO_LONDON}
      >
        Zhao and London (Transmembrane Tendency)
      </option>
    );
  }

  return (
    <ColorModeSelectLabel>
      Color by:
      <select
        value={plotSettingsStore.tableColorMode}
        onChange={handleColorModeChange}
      >
        {colorModeElements}
      </select>
    </ColorModeSelectLabel>
  );
});

const CompareModeSelectLabel = styled.label`
  select {
    margin-left: 0.65em;
    max-width: 150px;
    border-radius: 3px;
  }
`;

const CompareModeSelect = observer(() => {
  const { plotSettingsStore, configStore } = useStores();

  const handleCompareModeChange = (event) =>
    plotSettingsStore.setTableCompareMode(event.target.value);
  const handleCompareColorChange = (event) =>
    plotSettingsStore.setTableCompareColor(event.target.value);

  const colorOptions = Object.keys(snapGeneHighlightColors);
  const colorOptionElements = [];
  // Add the color-code option
  colorOptionElements.push(
    <option
      key={COMPARE_COLORS.COLOR_MODE_CODE}
      value={COMPARE_COLORS.COLOR_MODE_CODE}
    >
      {configStore.dnaOrAa === 'dna' ? '4' : '20'}-color code
    </option>
  );

  if (configStore.dnaOrAa === 'aa') {
    colorOptionElements.push(
      <option
        key={COMPARE_COLORS.COLOR_MODE_CLUSTAL}
        value={COMPARE_COLORS.COLOR_MODE_CLUSTAL}
      >
        ClustalX (Properties)
      </option>,
      <option
        key={COMPARE_COLORS.COLOR_MODE_ZAPPO}
        value={COMPARE_COLORS.COLOR_MODE_ZAPPO}
      >
        Zappo (Physiochemical Properties)
      </option>,
      <option
        key={COMPARE_COLORS.COLOR_MODE_ZHAO_LONDON}
        value={COMPARE_COLORS.COLOR_MODE_ZHAO_LONDON}
      >
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
  if (configStore.dnaOrAa === 'aa') {
    colorOptionElements.push(
      <option
        key={COMPARE_COLORS.COMPARE_COLOR_DOTS}
        value={COMPARE_COLORS.COMPARE_COLOR_DOTS}
      >
        Dots
      </option>
    );
  }

  const disabled =
    plotSettingsStore.tableColorMode === COLOR_MODES.COLOR_MODE_CODE;

  return (
    <CompareModeSelectLabel>
      Color bases that:
      <select
        disabled={disabled}
        value={plotSettingsStore.tableCompareMode}
        onChange={handleCompareModeChange}
      >
        <option value={COMPARE_MODES.COMPARE_MODE_MATCH}>Match</option>
        <option value={COMPARE_MODES.COMPARE_MODE_MISMATCH}>
          Don&apos;t Match
        </option>
      </select>
      <select
        disabled={disabled}
        value={plotSettingsStore.tableCompareColor}
        onChange={handleCompareColorChange}
      >
        {colorOptionElements}
      </select>
    </CompareModeSelectLabel>
  );
});

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

const TableOptions = observer(() => {
  const { dataStore } = useStores();

  const handleDownloadSelect = (option) => {
    if (option === 'Acknowledgements') {
      dataStore.downloadAcknowledgements();
    } else if (option === 'Aggregate Data') {
      dataStore.downloadAggCaseData();
    } else if (option === 'Accession IDs') {
      dataStore.downloadAccessionIds();
    }
  };

  return (
    <DataTableOptions>
      <ColorModeSelect />
      <CompareModeSelect />
      <Spacer />
      <DropdownButton
        text={'Download'}
        options={['Aggregate Data', 'Acknowledgements', 'Accession IDs']}
        onSelect={handleDownloadSelect}
      />
    </DataTableOptions>
  );
});
export default TableOptions;
