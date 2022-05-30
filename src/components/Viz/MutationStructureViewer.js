import React, { useState, useRef, useEffect } from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import './../../styles/litemol.min.css';

import LiteMol from 'litemol';
const Bootstrap = LiteMol.Bootstrap;
const Transformer = Bootstrap.Entity.Transformer;

import { reds } from '../../constants/colors';
import { hexToRgb } from '../../utils/color';
import {
  colorHeatmap,
  getMoleculeAssemblies,
  getMoleculeEntities,
  CreateMacromoleculeVisual,
} from '../LiteMol/litemolutils';
import { getProtein } from '../../utils/gene_protein';
import {
  COORDINATE_MODES,
  LITEMOL_STYLES,
  NORM_MODES,
} from '../../constants/defs.json';

import EmptyPlot from '../Common/EmptyPlot';
import DropdownButton from '../Buttons/DropdownButton';
import DownloadPymolScriptModal from '../Modals/DownloadPymolScriptModal';
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
} from './MutationStructureViewer.styles';

const NUM_COLORS = reds.length;
const DOWNLOAD_OPTIONS = {
  DOWNLOAD_DATA: 'Download Data',
  DOWNLOAD_PYMOL: 'Download PyMOL Script',
};

const MutationStructureViewer = observer(() => {
  const { configStore, dataStore, plotSettingsStore } = useStores();
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
    plotSettingsStore.setMutationStructureActiveAssembly(event.target.value);
  };
  const onChangeProteinStyle = (event) => {
    plotSettingsStore.setMutationStructureProteinStyle(event.target.value);
  };
  const onChangeNormMode = (event) => {
    plotSettingsStore.setMutationStructureNormMode(event.target.value);
  };

  const applyChanges = () => {
    plotSettingsStore.setReportStructurePdbId(state.pdbId);

    // Clear changed and error states
    setState({
      ...state,
      pdbIdChanged: false,
      validPdbId: true,
    });
  };

  const handleDownloadSelect = (downloadOption) => {
    if (downloadOption === DOWNLOAD_OPTIONS.DOWNLOAD_DATA) {
      // groupDataStore.downloadStructureMutationData();
      console.log('DOWNLOAD DATA');
    } else if (downloadOption === DOWNLOAD_OPTIONS.DOWNLOAD_PYMOL) {
      showDownloadPymolScriptModal();
    }
  };

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
        <EmptyPlot height={250}>
          <p>
            No analogous protein found for gene {configStore.selectedGene.name}.
            Please click on the "Filter Sequences" button at the top of the
            screen and switch to "Protein" under "Genomic Coordinates" to enter
            protein mode.
          </p>
        </EmptyPlot>
      );
    }
  }

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
        mutation.colorInd = Math.floor(
          (mutation[colorField] - 0.001) * NUM_COLORS
        );
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
    plugin.clear();

    const pdbId = plotSettingsStore.mutationStructurePdbId.toLowerCase();

    // good example: https://github.com/dsehnal/LiteMol/blob/master/src/Viewer/App/Examples.ts
    const modelAction = plugin
      .createTransform()
      .add(plugin.root, Transformer.Data.Download, {
        url: `https://www.ebi.ac.uk/pdbe/static/entry/${pdbId}_updated.cif`,
        type: 'String',
        id: pdbId,
      })
      .then(Transformer.Data.ParseCif, { id: pdbId }, { isBinding: true })
      .then(
        Transformer.Molecule.CreateFromMmCif,
        { blockIndex: 0 },
        { ref: 'molecule' }
      );

    plugin.applyTransform(modelAction).then(() => {
      let vizAction = plugin
        .createTransform()
        .add(
          'molecule',
          Transformer.Molecule.CreateModel,
          { modelIndex: 0 },
          { ref: 'model' }
        );

      // If an assembly exists, then display that instead
      // of the asymmetric unit
      const assemblies = getMoleculeAssemblies({ plugin });

      if (assemblies.length > 0) {
        // If no assembly is selected, then default to the first assembly
        if (useAssembly === '') {
          useAssembly = assemblies[0];
        }
        //useAssembly = 'asym';

        // If user decides to display the asymmetric unit,
        // then skip the assembly process
        if (useAssembly !== 'asym') {
          vizAction = vizAction.then(
            Transformer.Molecule.CreateAssembly,
            { name: assemblies[0] },
            { ref: 'assembly' }
          );
        }
      }

      const entities = getMoleculeEntities({ plugin });

      setState({
        ...state,
        assemblies,
        entities,
        activeAssembly: useAssembly,
      });

      vizAction = vizAction.then(CreateMacromoleculeVisual, {
        polymer: true,
        het: true,
        water: false,
        style: plotSettingsStore.mutationStructureProteinStyle,
      });

      plugin.applyTransform(vizAction).then(() => {
        applyHeatmap({
          ref:
            assemblies.length > 0 && useAssembly !== 'asym'
              ? 'assembly'
              : 'model',
          entities: entities,
        });
      });

      // Update store so other components have this info
      plotSettingsStore.setMutationStructureAssemblies(assemblies);
      plotSettingsStore.setMutationStructureActiveAssembly(useAssembly);
      plotSettingsStore.setMutationStructureEntities(entities);
    });

    modelAction;
  };

  useEffect(() => {
    if (!plugin) return;
    loadModel();
  }, [
    plugin,
    plotSettingsStore.mutationStructurePdbId,
    plotSettingsStore.mutationStructureProteinStyle,
  ]);

  useEffect(() => {
    if (!plugin) return;
    applyHeatmap({
      ref:
        state.assemblies.length > 0 && state.activeAssembly !== 'asym'
          ? 'assembly'
          : 'model',
      entities: state.entities,
    });
  }, [plotSettingsStore.mutationStructureNormMode]);

  const assemblyOptionItems = [];
  state.assemblies.forEach((assembly) => {
    assemblyOptionItems.push(
      <option key={`assembly-option-${assembly}`} value={assembly}>
        {assembly}
      </option>
    );
  });

  return (
    <MutationStructureViewerContainer>
      <DownloadPymolScriptModal
        isOpen={state.downloadPymolScriptModalOpen}
        onRequestClose={hideDownloadPymolScriptModal}
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
      <LiteMolContainer>
        <LiteMolPlugin
          height={500}
          plugin={plugin}
          ref={pluginRef}
          setPlugin={setPlugin}
        />
      </LiteMolContainer>
    </MutationStructureViewerContainer>
  );
});

export default MutationStructureViewer;
