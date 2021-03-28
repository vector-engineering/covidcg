import React, { useState, useEffect, useRef } from 'react';
import './../../styles/litemol.min.css';

import LiteMol from 'litemol';

const Rx = LiteMol.Core.Rx;
const Bootstrap = LiteMol.Bootstrap;
const Transformer = Bootstrap.Entity.Transformer;

const LiteMolPlugin = React.forwardRef(({ plugin, setPlugin }, ref) => {
  useEffect(() => {
    const _plugin = LiteMol.Plugin.create({
      target: ref.current,
      layoutState: {
        hideControls: true,
        collapsedControlsLayout:
          Bootstrap.Components.CollapsedControlsLayout.Landscape,
      },
      viewportBackground: '#F1F1F1',
    });
    setPlugin(_plugin);
  }, []);

  return (
    <div
      ref={ref}
      style={{
        position: 'absolute',
        top: 24,
        right: 0,
        left: 0,
        bottom: 0,
      }}
    />
  );
});

const LiteMolTab = () => {
  const [plugin, setPlugin] = useState(null);
  const pluginRef = useRef(null);

  function createTheme(colorDef) {
    const model = Bootstrap.Tree.Selection.byRef(pluginRef)
      .subtree()
      .ofType(Bootstrap.Entity.Molecule.Visual);
    console.log(model);
    const mapper = new ColorMapper();
    mapper.addColor(colorDef.base);
    const map = new Uint8Array(model.data.atoms.count);

    for (const e of colorDef.entries) {
      const query = Q.sequence(
        e.entity_id.toString(),
        e.struct_asym_id,
        { seqNumber: e.start_residue_number },
        { seqNumber: e.end_residue_number }
      ).compile();
      const colorIndex = mapper.addColor(e.color);
      for (const f of query(model.queryContext).fragments) {
        for (const a of f.atomIndices) {
          map[a] = colorIndex;
        }
      }
    }

    const fallbackColor = { r: 0.6, g: 0.6, b: 0.6 };
    const selectionColor = { r: 0, g: 0, b: 1 };
    const highlightColor = { r: 1, g: 0, b: 1 };

    const colors = Core.Utils.FastMap.create();
    colors.set('Uniform', fallbackColor);
    colors.set('Selection', selectionColor);
    colors.set('Highlight', highlightColor);

    const mapping = Visualization.Theme.createColorMapMapping(
      (i) => map[i],
      mapper.colorMap,
      fallbackColor
    );
    // make the theme "sticky" so that it persist "ResetScene" command.
    return Visualization.Theme.createMapping(mapping, {
      colors,
      isSticky: true,
    });
  }

  function applyTheme(theme) {
    const visuals = plugin.selectEntities(
      Bootstrap.Tree.Selection.byRef(pluginRef)
        .subtree()
        .ofType(Bootstrap.Entity.Molecule.Visual)
    );
    for (const v of visuals) {
      plugin.command(Bootstrap.Command.Visual.UpdateBasicTheme, {
        visual: v,
        theme,
      });
    }
  }

  const colorSequences = () => {
    const coloring = {
      base: { r: 255, g: 255, b: 255 },
      entries: [
        {
          entity_id: '1',
          struct_asym_id: 'A',
          start_residue_number: 319,
          end_residue_number: 541,
          color: { r: 255, g: 128, b: 64 },
        },
        {
          entity_id: '1',
          struct_asym_id: 'A',
          start_residue_number: 14,
          end_residue_number: 305,
          color: { r: 64, g: 128, b: 255 },
        },
      ],
    };

    const theme = createTheme(coloring);
    // instead of "polymer-visual", "model" or any valid ref can be used: all "child" visuals will be colored.
    applyTheme(theme);
  };

  useEffect(() => {
    if (!plugin) {
      return;
    }
    let action = plugin.createTransform();
    const pdbId = '6zgg';

    action
      .add(plugin.root, Transformer.Data.Download, {
        url: `https://www.ebi.ac.uk/pdbe/static/entry/${pdbId}_updated.cif`,
        type: 'String',
        id: pdbId,
      })
      .then(Transformer.Data.ParseCif, { id: pdbId }, { isBinding: true })
      .then(Transformer.Molecule.CreateFromMmCif, { blockIndex: 0 })
      .then(Transformer.Molecule.CreateModel, { modelIndex: 0 })
      .then(Transformer.Molecule.CreateMacromoleculeVisual, {
        polymer: true,
        het: true,
        water: false,
      });

    plugin.applyTransform(action);
  }, [plugin]);

  return (
    <>
      <LiteMolPlugin plugin={plugin} ref={pluginRef} setPlugin={setPlugin} />
      <button onClick={colorSequences}>color sequences</button>
    </>
  );
};

export default LiteMolTab;
