import React, { useState, useRef, useEffect } from 'react';
import LiteMol from 'litemol';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import './../../styles/litemol.min.css';

import {
  colorHeatmap,
  getMoleculeAssemblies,
  CreateMacromoleculeVisual,
} from '../LiteMol/litemolutils';
import { reds } from '../../constants/colors';
import { LITEMOL_STYLES } from '../../constants/defs';
import { hexToRgb } from '../../utils/color';
import { getAllProteins } from '../../utils/gene_protein';
import defaultStructures from '../../../static_data/default_structures.json';

import EmptyPlot from '../Common/EmptyPlot';
import DownloadPymolScriptModal from '../Modals/DownloadPymolScriptModal';
import LiteMolPlugin from '../LiteMol/LiteMolPlugin';
import {
  StructuralViewerContainer,
  LiteMolContainer,
  StructuralViewerHeader,
  OptionSelectContainer,
  OptionInputContainer,
  ConfirmButton,
  InvalidText,
} from './GroupStructuralViewer.styles';

const proteins = getAllProteins();

const Bootstrap = LiteMol.Bootstrap;
const Transformer = Bootstrap.Entity.Transformer;

const numColors = reds.length;

const StructuralViewer = observer(() => {
  const { groupDataStore, plotSettingsStore } = useStores();
  const [plugin, setPlugin] = useState(null);
  const [state, setState] = useState({
    downloadPymolScriptModalOpen: false,
    activeProtein: plotSettingsStore.reportStructureActiveProtein,
    activeProteinChanged: false,
    pdbId: plotSettingsStore.reportStructurePdbId,
    validPdbId: true,
    pdbIdChanged: false,
    assemblies: [],
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

  const onChangeStructureActiveGroup = (event) => {
    plotSettingsStore.setReportStructureActiveGroup(event.target.value);
  };

  const onChangeReportStructureActiveProtein = (event) => {
    const newProtein = event.target.value;

    if (newProtein === state.activeProtein) {
      return;
    }

    setState({
      ...state,
      activeProtein: newProtein,
      activeProteinChanged:
        newProtein != plotSettingsStore.reportStructureActiveProtein,
      // Clear PDB ID when protein changes
      pdbId: Object.keys(defaultStructures).includes(newProtein)
        ? defaultStructures[newProtein]
        : '',
      validPdbId: Object.keys(defaultStructures).includes(newProtein)
        ? true
        : false,
      pdbIdChanged: true,
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
      pdbIdChanged: newPdbId != plotSettingsStore.reportStructurePdbId,
    });
  };

  const onChangeActiveAssembly = (event) => {
    loadModel({ useAssembly: event.target.value });
  };

  const onChangeProteinStyle = (event) => {
    plotSettingsStore.setReportStructureProteinStyle(event.target.value);
  };

  const applyChanges = () => {
    plotSettingsStore.setReportStructureActiveProtein(state.activeProtein);
    plotSettingsStore.setReportStructurePdbId(state.pdbId);

    // Clear changed and error states
    setState({
      ...state,
      activeProteinChanged: false,
      pdbIdChanged: false,
      validPdbId: true,
    });
  };

  const downloadData = () => {
    groupDataStore.downloadStructureMutationData();
  };

  const applyHeatmap = ({ ref }) => {
    const mutations = groupDataStore.groupMutationFrequency[
      groupDataStore.activeGroupType
    ]['protein_aa']['0']
      .filter(
        (groupMutation) =>
          groupMutation.name === plotSettingsStore.reportStructureActiveGroup &&
          groupMutation.protein ===
            plotSettingsStore.reportStructureActiveProtein
      )
      // Convert fractional frequencies to colors
      .slice()
      .map((mut) => {
        mut.colorInd = Math.floor((mut.fraction - 0.001) * numColors);
        return mut;
      })
      .sort((a, b) => a.pos - b.pos);

    // For each color level, build a list of residues
    const heatmapEntries = [];
    reds.forEach((color, i) => {
      let indices = mutations
        .filter((mut) => mut.colorInd == i)
        .map((mut) => mut.pos);
      if (indices.length === 0) return;

      heatmapEntries.push({
        indices,
        color: hexToRgb(color.substr(1)),
      });
    });

    colorHeatmap({ plugin, entries: heatmapEntries, ref });
  };

  const loadModel = ({ useAssembly } = { useAssembly: '' }) => {
    if (!plugin) {
      return;
    }
    plugin.clear();

    const pdbId = plotSettingsStore.reportStructurePdbId.toLowerCase();

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

        // If user decides to display the asymmetric unit,
        // then skip the assembly process
        if (useAssembly !== 'asym') {
          vizAction = vizAction.then(
            Transformer.Molecule.CreateAssembly,
            { name: assemblies[0] },
            { ref: 'assembly' }
          );
        }

        // TODO: remove the original model from the tree?
        setState({
          ...state,
          assemblies,
          activeAssembly: useAssembly,
        });
      }

      vizAction = vizAction.then(CreateMacromoleculeVisual, {
        polymer: true,
        het: true,
        water: false,
        style: plotSettingsStore.reportStructureProteinStyle,
      });

      plugin.applyTransform(vizAction).then(() => {
        applyHeatmap({
          ref:
            assemblies.length > 0 && useAssembly !== 'asym'
              ? 'assembly'
              : 'model',
        });
      });
    });

    modelAction;
  };

  useEffect(() => {
    loadModel();
  }, [
    plugin,
    plotSettingsStore.reportStructurePdbId,
    plotSettingsStore.reportStructureActiveProtein,
    plotSettingsStore.reportStructureProteinStyle,
  ]);

  useEffect(() => {
    if (!plugin) return;
    applyHeatmap({
      ref:
        state.assemblies.length > 0 && state.activeAssembly !== 'asym'
          ? 'assembly'
          : 'model',
    });
  }, [plotSettingsStore.reportStructureActiveGroup]);

  const proteinOptionItems = [];
  proteins.forEach((protein) => {
    proteinOptionItems.push(
      <option key={`structure-protein-${protein.name}`} value={protein.name}>
        {protein.name}
      </option>
    );
  });

  const groupOptionItems = [];
  groupDataStore.selectedGroups.forEach((group) => {
    groupOptionItems.push(
      <option key={`structure-active-group-${group}`} value={group}>
        {group}
      </option>
    );
  });

  if (groupDataStore.selectedGroups.length === 0) {
    return (
      <EmptyPlot height={250}>
        <p>No {groupDataStore.getActiveGroupTypePrettyName()}s selected</p>
      </EmptyPlot>
    );
  }

  const assemblyOptionItems = [];
  state.assemblies.forEach((assembly) => {
    assemblyOptionItems.push(
      <option key={`assembly-option-${assembly}`} value={assembly}>
        {assembly}
      </option>
    );
  });

  return (
    <StructuralViewerContainer>
      <DownloadPymolScriptModal
        isOpen={state.downloadPymolScriptModalOpen}
        onRequestClose={hideDownloadPymolScriptModal}
      />
      <StructuralViewerHeader>
        <OptionSelectContainer>
          <label>
            Displaying mutations for{' '}
            {groupDataStore.getActiveGroupTypePrettyName()}
            <select
              value={plotSettingsStore.reportStructureActiveGroup}
              onChange={onChangeStructureActiveGroup}
            >
              {groupOptionItems}
            </select>
          </label>
        </OptionSelectContainer>
        <div className="spacer"></div>
        <ConfirmButton onClick={downloadData}>Download Data</ConfirmButton>
      </StructuralViewerHeader>
      <StructuralViewerHeader>
        <OptionInputContainer>
          <label>
            PDB ID
            <input type="text" value={state.pdbId} onChange={onChangePdbId} />
          </label>
          {!state.validPdbId && <InvalidText>Invalid PDB ID</InvalidText>}
          {(state.pdbIdChanged || state.activeProteinChanged) && (
            <ConfirmButton disabled={!state.validPdbId} onClick={applyChanges}>
              Apply
            </ConfirmButton>
          )}
        </OptionInputContainer>
        <OptionSelectContainer>
          <label>
            Protein
            <select
              value={state.activeProtein}
              onChange={onChangeReportStructureActiveProtein}
            >
              {proteinOptionItems}
            </select>
          </label>
        </OptionSelectContainer>
        <div className="spacer"></div>
        <ConfirmButton onClick={showDownloadPymolScriptModal}>
          Download PyMOL Script
        </ConfirmButton>
      </StructuralViewerHeader>
      <StructuralViewerHeader>
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
              value={plotSettingsStore.reportStructureProteinStyle}
              onChange={onChangeProteinStyle}
            >
              <option value={LITEMOL_STYLES.SURFACE}>Surface</option>
              <option value={LITEMOL_STYLES.CARTOON}>Cartoon</option>
            </select>
          </label>
        </OptionSelectContainer>
      </StructuralViewerHeader>
      <LiteMolContainer>
        <LiteMolPlugin
          height={500}
          plugin={plugin}
          ref={pluginRef}
          setPlugin={setPlugin}
        />
      </LiteMolContainer>
    </StructuralViewerContainer>
  );
});

export default StructuralViewer;
