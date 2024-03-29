import React, { useState, useRef, useEffect } from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import './../../styles/litemol.min.css';

import { reds } from '../../constants/colors';
import { hexToRgb } from '../../utils/color';
import { LoadLitemolModel, colorHeatmap } from '../LiteMol/litemolutils';
import { getProtein } from '../../utils/gene_protein';
import {
  DNA_OR_AA,
  COORDINATE_MODES,
  LITEMOL_STYLES,
  NORM_MODES,
} from '../../constants/defs.json';

import EmptyPlot from '../Common/EmptyPlot';
import DropdownButton from '../Buttons/DropdownButton';
import DownloadPymolScriptModal from '../Modals/DownloadPymolScriptModal';
import StructureEntities from '../LiteMol/StructureEntities';
import LiteMolPlugin from '../LiteMol/LiteMolPlugin';

import {
  PlotOptions,
  OptionSelectContainer,
  OptionInputContainer,
} from './Plot.styles';
import {
  MutationStructureViewerContainer,
  LiteMolContainer,
  InvalidText,
  ConfirmButton,
  HighlightedMutations,
} from './MutationStructureViewer.styles';

const NUM_COLORS = reds.length;
const DOWNLOAD_OPTIONS = {
  DOWNLOAD_DATA: 'Download Data',
  DOWNLOAD_PYMOL: 'Download PyMOL Script',
};

const MutationStructureViewer = observer(() => {
  const { configStore, dataStore, mutationDataStore, plotSettingsStore } =
    useStores();
  const [plugin, setPlugin] = useState(null);
  const [state, setState] = useState({
    downloadPymolScriptModalOpen: false,
    pdbId: plotSettingsStore.reportStructurePdbId,
    validPdbId: true,
    pdbIdChanged: false,
    assemblies: [],
    entities: [],
    activeAssembly: '',
  });
  const pluginRef = useRef(null);

  // Flag to not run any hooks if we're going to render an empty plot
  // i.e., show the hint text when in NT mode
  const skipHooks =
    configStore.dnaOrAa !== DNA_OR_AA.AA ||
    (configStore.coordinateMode !== COORDINATE_MODES.COORD_GENE &&
      configStore.coordinateMode !== COORDINATE_MODES.COORD_PROTEIN);

  const showDownloadPymolScriptModal = () => {
    setState({
      ...state,
      downloadPymolScriptModalOpen: true,
    });
  };
  const hideDownloadPymolScriptModal = () => {
    setState({
      ...state,
      downloadPymolScriptModalOpen: false,
    });
  };

  const onChangePdbId = (event) => {
    const newPdbId = event.target.value;
    let validPdbId = true;

    if (newPdbId.length !== 4) {
      validPdbId = false;
    }

    setState({
      ...state,
      pdbId: newPdbId,
      validPdbId,
      pdbIdChanged: newPdbId != plotSettingsStore.mutationStructurePdbId,
    });
  };

  const onChangeActiveAssembly = (event) => {
    loadModel({ useAssembly: event.target.value });
    plotSettingsStore.applyPendingChanges({
      mutationStructureActiveAssembly: event.target.value,
    });
  };
  const onChangeProteinStyle = (event) => {
    plotSettingsStore.applyPendingChanges({
      mutationStructureProteinStyle: event.target.value,
    });
  };
  const onChangeNormMode = (event) => {
    plotSettingsStore.applyPendingChanges({
      mutationStructureNormMode: event.target.value,
    });
  };

  const applyChanges = () => {
    plotSettingsStore.applyPendingChanges({
      mutationStructurePdbId: state.pdbId,
    });

    // Clear changed and error states
    setState({
      ...state,
      pdbIdChanged: false,
      validPdbId: true,
    });
  };

  const handleDownloadSelect = (downloadOption) => {
    if (downloadOption === DOWNLOAD_OPTIONS.DOWNLOAD_DATA) {
      dataStore.downloadMutationFrequencies();
    } else if (downloadOption === DOWNLOAD_OPTIONS.DOWNLOAD_PYMOL) {
      showDownloadPymolScriptModal();
    }
  };
  const handlePymolScriptDownload = (opts) => {
    dataStore.downloadMutationStructurePymolScript(opts);
  };

  const onChangeEntities = (entities) => {
    setState({
      ...state,
      entities,
    });
    applyHeatmap({
      ref:
        state.assemblies.length > 0 && state.activeAssembly !== 'asym'
          ? 'assembly'
          : 'model',
      entities,
    });
    plotSettingsStore.applyPendingChanges({
      mutationStructureEntities: entities,
    });
  };

  const applyHeatmap = ({ ref, entities }) => {
    // Deep copy data from store
    const mutations = JSON.parse(JSON.stringify(dataStore.groupCounts))
      // Convert fractional frequencies to colors
      .map((mutation) => {
        let colorField;
        if (
          plotSettingsStore.mutationStructureNormMode ===
          NORM_MODES.NORM_PERCENTAGES
        ) {
          colorField = 'percent';
        } else if (
          plotSettingsStore.mutationStructureNormMode ===
          NORM_MODES.NORM_COVERAGE_ADJUSTED
        ) {
          colorField = 'partial_adjusted';
        }
        mutation.colorInd = Math.floor(mutation[colorField] * (NUM_COLORS - 1));
        return mutation;
      })
      // Sort by position
      .sort((a, b) => a.pos - b.pos);

    // console.log(mutations);

    // For each color level, build a list of residues
    const heatmapEntries = [];
    reds.forEach((color, i) => {
      let indices = mutations
        .filter((mut) => mut.colorInd == i)
        .map((mut) => mut.pos);
      if (indices.length === 0) return;

      heatmapEntries.push({
        indices,
        color: hexToRgb(color.substring(1)),
      });
    });

    // Add selected mutations
    const selectedMutationPositions = [];
    configStore.selectedGroups.forEach(({ group }) => {
      const mutation = mutationDataStore.mutationStrToMutation(
        configStore.dnaOrAa,
        configStore.coordinateMode,
        group
      );

      if (mutation !== undefined) {
        selectedMutationPositions.push(mutation.pos);
      }
    });

    if (selectedMutationPositions.length > 0) {
      heatmapEntries.push({
        indices: selectedMutationPositions,
        color: hexToRgb('0000ff'),
      });
    }

    const ignoreChains = [];
    entities.forEach((entity) => {
      if (!entity.checked) {
        ignoreChains.push(...entity.chains);
      }
    });

    colorHeatmap({ plugin, entries: heatmapEntries, ref, ignoreChains });
  };

  const loadModel = ({ useAssembly } = { useAssembly: '' }) => {
    if (!plugin) {
      return;
    }

    LoadLitemolModel({
      plugin,
      pdbId: plotSettingsStore.mutationStructurePdbId,
      proteinStyle: plotSettingsStore.mutationStructureProteinStyle,
      useAssembly,
      onLoad: ({ vizAction, assemblies, entities, activeAssembly }) => {
        setState({
          ...state,
          assemblies,
          entities,
          activeAssembly,
        });

        plugin.applyTransform(vizAction).then(() => {
          applyHeatmap({
            ref:
              assemblies.length > 0 && activeAssembly !== 'asym'
                ? 'assembly'
                : 'model',
            entities: entities,
          });
        });

        // Update store so other components have this info
        plotSettingsStore.applyPendingChanges({
          mutationStructureAssemblies: assemblies,
          mutationStructureActiveAssembly: activeAssembly,
          mutationStructureEntities: entities,
        });
      },
    });
  };

  useEffect(() => {
    if (skipHooks) return;
    if (!plugin) return;
    loadModel();
  }, [
    plugin,
    plotSettingsStore.mutationStructurePdbId,
    plotSettingsStore.mutationStructureProteinStyle,
  ]);

  useEffect(() => {
    if (skipHooks) return;
    if (!plugin) return;
    applyHeatmap({
      ref:
        state.assemblies.length > 0 && state.activeAssembly !== 'asym'
          ? 'assembly'
          : 'model',
      entities: state.entities,
    });
  }, [plotSettingsStore.mutationStructureNormMode, configStore.selectedGroups]);

  const assemblyOptionItems = [];
  state.assemblies.forEach((assembly) => {
    assemblyOptionItems.push(
      <option key={`assembly-option-${assembly}`} value={assembly}>
        {assembly}
      </option>
    );
  });

  // Only valid in AA mode, or gene/protein mode
  if (configStore.dnaOrAa !== DNA_OR_AA.AA) {
    return (
      <EmptyPlot height={150}>
        <p>
          Structure view only available in AA mode. Please select the
          &quot;AA&quot; option under &quot;Mutation Format&quot; in the top bar
          or in the &quot;Filter Sequences&quot; dialog.
        </p>
      </EmptyPlot>
    );
  } else if (
    configStore.coordinateMode !== COORDINATE_MODES.COORD_GENE &&
    configStore.coordinateMode !== COORDINATE_MODES.COORD_PROTEIN
  ) {
    return (
      <EmptyPlot height={150}>
        <p>
          Structure view requires gene or protein coordinates. Please select the
          &quot;Gene&quot; or &quot;Protein&quot; option under &quot;Coordinate
          Mode&quot; in the top bar or in the &quot;Filter Sequences&quot;
          dialog.
        </p>
      </EmptyPlot>
    );
  }

  // Get the protein to show
  let activeProtein;
  if (configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
    activeProtein = configStore.selectedProtein;
  }
  // If we're in gene mode, we need to find the analogous protein
  // If it doesn't exist (i.e., gene is an ORF that codes for multiple proteins),
  // then show a message prompting user to switch to protein mode
  else if (configStore.coordinateMode === COORDINATE_MODES.COORD_GENE) {
    activeProtein = getProtein(
      configStore.selectedGene.name,
      configStore.selectedReference
    );
    // If we fell back to another protein, then the analogous protein doesn't exist
    if (activeProtein.name !== configStore.selectedGene.name) {
      return (
        <EmptyPlot height={150}>
          <p>
            No analogous protein found for gene {configStore.selectedGene.name}.
            Please click on the &quot;Filter Sequences&quot; button at the top
            of the screen and switch to &quot;Protein&quot; under &quot;Genomic
            Coordinates&quot; to enter protein mode.
          </p>
        </EmptyPlot>
      );
    }
  }

  // Selected mutations, sort by position
  const selectedMutationItems = [];
  JSON.parse(JSON.stringify(configStore.selectedGroups))
    .map(({ group }) => {
      const mutation = mutationDataStore.mutationStrToMutation(
        configStore.dnaOrAa,
        configStore.coordinateMode,
        group
      );
      return mutation;
    })
    .filter((mutation) => mutation !== undefined)
    .sort((a, b) => {
      return a.pos - b.pos;
    })
    .forEach((mutation) => {
      selectedMutationItems.push(
        <li
          key={`mutation-structure-selected-mutation-${mutation.mutation_str}`}
        >
          {mutation.name}
        </li>
      );
    });

  return (
    <MutationStructureViewerContainer>
      <DownloadPymolScriptModal
        isOpen={state.downloadPymolScriptModalOpen}
        onRequestClose={hideDownloadPymolScriptModal}
        onConfirm={handlePymolScriptDownload}
      />
      <PlotOptions style={{ marginBottom: '2px' }}>
        <OptionInputContainer>
          <label>
            PDB ID
            <input type="text" value={state.pdbId} onChange={onChangePdbId} />
          </label>
          {!state.validPdbId && <InvalidText>Invalid PDB ID</InvalidText>}
        </OptionInputContainer>
        {state.pdbIdChanged && (
          <ConfirmButton disabled={!state.validPdbId} onClick={applyChanges}>
            Apply
          </ConfirmButton>
        )}
        <OptionSelectContainer>
          <label>
            Show:
            <select
              value={plotSettingsStore.mutationStructureNormMode}
              onChange={onChangeNormMode}
            >
              <option value={NORM_MODES.NORM_PERCENTAGES}>Percents</option>
              <option value={NORM_MODES.NORM_COVERAGE_ADJUSTED}>
                Percents (coverage adjusted)
              </option>
            </select>
          </label>
        </OptionSelectContainer>
        <div className="spacer" />
        <DropdownButton
          text={'Download'}
          options={[
            DOWNLOAD_OPTIONS.DOWNLOAD_DATA,
            DOWNLOAD_OPTIONS.DOWNLOAD_PYMOL,
          ]}
          onSelect={handleDownloadSelect}
        />
      </PlotOptions>
      <PlotOptions>
        <OptionSelectContainer>
          <label>
            Assembly
            <select
              value={state.activeAssembly}
              onChange={onChangeActiveAssembly}
            >
              <option value="asym">Asymmetric Unit</option>
              {assemblyOptionItems}
            </select>
          </label>
        </OptionSelectContainer>
        <OptionSelectContainer>
          <label>
            Protein Style
            <select
              value={plotSettingsStore.mutationStructureProteinStyle}
              onChange={onChangeProteinStyle}
            >
              <option value={LITEMOL_STYLES.SURFACE}>Surface</option>
              <option value={LITEMOL_STYLES.CARTOON}>Cartoon</option>
            </select>
          </label>
        </OptionSelectContainer>
      </PlotOptions>
      <StructureEntities
        entities={state.entities}
        onChangeEntities={onChangeEntities}
      />
      <LiteMolContainer>
        <LiteMolPlugin
          height={500}
          plugin={plugin}
          ref={pluginRef}
          setPlugin={setPlugin}
        />
      </LiteMolContainer>
      {selectedMutationItems.length > 0 && (
        <HighlightedMutations>
          Highlighted mutations:
          <ul>{selectedMutationItems}</ul>
        </HighlightedMutations>
      )}
    </MutationStructureViewerContainer>
  );
});

export default MutationStructureViewer;
