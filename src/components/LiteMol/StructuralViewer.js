import React, { useState, useRef, useEffect } from 'react';
import LiteMol from 'litemol';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import './../../styles/litemol.min.css';

import { colorHeatmap } from './litemolutils';
import { reds } from '../../constants/colors';
import { hexToRgb } from '../../utils/color';

import LiteMolPlugin from './LiteMolPlugin';

import {
  StructuralViewerContainer,
  LiteMolContainer,
  StructuralViewerHeader,
  OptionSelectContainer,
  OptionInputContainer,
} from './StructuralViewer.styles';

const Core = LiteMol.Core;
const Bootstrap = LiteMol.Bootstrap;
const Transformer = Bootstrap.Entity.Transformer;
const Visualization = LiteMol.Visualization;

const numColors = reds.length;

const StructuralViewer = observer(() => {
  const { groupDataStore } = useStores();
  const [plugin, setPlugin] = useState(null);
  const pluginRef = useRef(null);

  const onChangeStructureActiveGroup = (event) => {
    groupDataStore.updateStructureActiveGroup(event.target.value);
  };

  const applyHeatmap = () => {
    const snvs = groupDataStore.groupSnvFrequency[
      groupDataStore.activeGroupType
    ][groupDataStore.groupSnvType]
      .filter(
        (groupSnv) =>
          groupSnv.name === groupDataStore.structureActiveGroup &&
          groupSnv.gene === 'S'
      )
      // Convert fractional frequencies to colors
      .slice()
      .map((snv) => {
        snv.colorInd = Math.floor((snv.fraction - 0.001) * numColors);
        return snv;
      });

    // For each color level, build a list of residues
    const heatmapEntries = [];
    reds.forEach((color, i) => {
      let indices = snvs
        .filter((snv) => snv.colorInd == i)
        .map((snv) => snv.pos);
      if (indices.length === 0) return;

      heatmapEntries.push({
        indices,
        color: hexToRgb(color.substr(1)),
      });
    });

    colorHeatmap({ plugin, entries: heatmapEntries });
  };

  useEffect(() => {
    if (!plugin) {
      return;
    }
    let action = plugin.createTransform();
    const pdbId = '6zgg';

    let modelAction = action
      .add(plugin.root, Transformer.Data.Download, {
        url: `https://www.ebi.ac.uk/pdbe/static/entry/${pdbId}_updated.cif`,
        type: 'String',
        id: pdbId,
      })
      .then(Transformer.Data.ParseCif, { id: pdbId }, { isBinding: true })
      .then(Transformer.Molecule.CreateFromMmCif, { blockIndex: 0 })
      .then(
        Transformer.Molecule.CreateModel,
        { modelIndex: 0 },
        { ref: 'model' }
      );
    // .then(Transformer.Molecule.CreateMacromoleculeVisual, {
    //   polymer: true,
    //   het: true,
    //   water: false,
    // });

    let selectionColors = Bootstrap.Immutable.Map()
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

    modelAction
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

    plugin.applyTransform(modelAction).then(() => {
      applyHeatmap();
    });
  }, [plugin]);

  useEffect(() => {
    if (!plugin) return;
    applyHeatmap();
  }, [groupDataStore.structureActiveGroup]);

  const groupOptionItems = [];
  groupDataStore.selectedGroups.forEach((group) => {
    groupOptionItems.push(
      <option key={`structure-active-group-${group}`} value={group}>
        {group}
      </option>
    );
  });

  return (
    <StructuralViewerContainer>
      <StructuralViewerHeader>
        <OptionSelectContainer>
          <label>
            <select
              value={groupDataStore.structureActiveGroup}
              onChange={onChangeStructureActiveGroup}
            >
              {groupOptionItems}
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
