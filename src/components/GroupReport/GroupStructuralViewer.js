import React, { useState, useRef, useEffect } from 'react';
import LiteMol from 'litemol';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import './../../styles/litemol.min.css';

import { colorHeatmap, getMoleculeAssemblies } from '../LiteMol/litemolutils';
import { reds } from '../../constants/colors';
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

const Core = LiteMol.Core;
const Bootstrap = LiteMol.Bootstrap;
const Transformer = Bootstrap.Entity.Transformer;
const Visualization = LiteMol.Visualization;

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
    activeAssembly: null,
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
    ]['protein_aa']
      .filter(
        (groupMutation) =>
          groupMutation.name === plotSettingsStore.reportStructureActiveGroup &&
          groupMutation.protein === plotSettingsStore.reportStructureActiveProtein
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

  const loadModel = () => {
    if (!plugin) {
      return;
    }
    plugin.clear();

    const pdbId = plotSettingsStore.reportStructurePdbId.toLowerCase();

    const selectionColors = Bootstrap.Immutable.Map()
      .set('Uniform', Visualization.Color.fromHex(0xaaaaaa))
      .set('Selection', Visualization.Theme.Default.SelectionColor)
      .set('Highlight', Visualization.Theme.Default.HighlightColor);
    const _style = {
      type: 'Surface',
      params: {
        probeRadius: 0,
        density: 1.25,
        smoothing: 3,
        isWireframe: false,
      },
      theme: {
        template: Bootstrap.Visualization.Molecule.Default.UniformThemeTemplate,
        colors: selectionColors,
        transparency: { alpha: 1.0 },
      },
    };

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
      // TODO: use the assembly information to allow the user
      // to select the chosen assembly, or the asymmetric unit
      const assemblies = getMoleculeAssemblies({ plugin });
      if (assemblies.length > 0) {
        vizAction = vizAction.then(
          Transformer.Molecule.CreateAssembly,
          { name: assemblies[0] },
          { ref: 'assembly' }
        );

        // TODO: remove the original model from the tree?

        setState({
          ...state,
          assemblies,
          activeAssembly: assemblies[0],
        });
      }

      vizAction = vizAction
        .then(
          Transformer.Molecule.CreateSelectionFromQuery,
          {
            query: Core.Structure.Query.nonHetPolymer(),
            name: 'Polymer',
            silent: true,
          },
          {}
        )
        .then(
          Transformer.Molecule.CreateVisual,
          { style: _style },
          { ref: 'Polymer' }
        );
      // .then(Transformer.Molecule.CreateMacromoleculeVisual, {
      //   polymer: true,
      //   het: true,
      //   water: false,
      // });

      plugin.applyTransform(vizAction).then(() => {
        applyHeatmap({ ref: assemblies.length > 0 ? 'assembly' : 'model' });
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
  ]);

  useEffect(() => {
    if (!plugin) return;
    applyHeatmap({ ref: state.assemblies.length > 0 ? 'assembly' : 'model' });
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

  return (
    <StructuralViewerContainer>
      <DownloadPymolScriptModal
        isOpen={state.downloadPymolScriptModalOpen}
        onRequestClose={hideDownloadPymolScriptModal}
      />
      <StructuralViewerHeader>
        <OptionSelectContainer>
          <label>
            Displaying mutations for {groupDataStore.getActiveGroupTypePrettyName()}
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
