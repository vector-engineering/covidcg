import React, { useState, useEffect, useRef } from 'react';
import './../../styles/litemol.min.css';

import LiteMol from 'litemol';

const Core = LiteMol.Core;
const Visualization = LiteMol.Visualization;
const Rx = LiteMol.Core.Rx;
const Bootstrap = LiteMol.Bootstrap;
const Transformer = Bootstrap.Entity.Transformer;

class ColorMapper {
  uniqueColors = [];
  map = Core.Utils.FastMap.create();

  get colorMap() {
    const map = Core.Utils.FastMap.create();
    this.uniqueColors.forEach((c, i) => map.set(i, c));
    return map;
  }

  addColor(color) {
    const id = `${color.r}-${color.g}-${color.b}`;
    if (this.map.has(id)) {
      return this.map.get(id);
    }
    const index = this.uniqueColors.length;
    this.uniqueColors.push(
      Visualization.Color.fromRgb(color.r, color.g, color.b)
    );
    this.map.set(id, index);
    return index;
  }
}

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
    let model = plugin.selectEntities('model');
    //   .subtree()
    //   .ofType(Bootstrap.Entity.Molecule.Visual);
    console.log(model);
    console.log(model.data);
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

  // todo: try this method
  const show_pdb = (id, regions, chain, entity) => {
    //this wont work cause its jquery
    $('#pdb_link').attr('href', 'https://www.ebi.ac.uk/pdbe/entry/pdb/' + id);
    var Bootstrap = LiteMol.Bootstrap;
    var Transformer = Bootstrap.Entity.Transformer;
    var Tree = Bootstrap.Tree;
    var Transform = Tree.Transform;
    LiteMol.Bootstrap.Command.Tree.RemoveNode.dispatch(
      plugin.context,
      plugin.context.tree.root
    );
    var action = Transform.build()
      .add(plugin.context.tree.root, Transformer.Data.Download, {
        url: 'https://www.ebi.ac.uk/pdbe/static/entry/' + id + '_updated.cif',
        type: 'String',
        id: id,
      })
      .then(Transformer.Data.ParseCif, { id: id }, { isBinding: true })
      .then(
        Transformer.Molecule.CreateFromMmCif,
        { blockIndex: 0 },
        { isBinding: true }
      )
      .then(
        Transformer.Molecule.CreateModel,
        { modelIndex: 0 },
        { isBinding: false, ref: 'model' }
      )
      .then(Transformer.Molecule.CreateMacromoleculeVisual, {
        polymer: true,
        polymerRef: 'polymer-visual',
        het: true,
        water: false,
      });
    applyTransforms(action).then(function (result) {
      var model = plugin.selectEntities('model')[0];
      if (!model) return;
      var coloring = {
        base: { r: 210, g: 180, b: 140 },
        entries: [],
      };
      $.each(regions.split(','), function (n, elem) {
        coloring['entries'].push({
          entity_id: entity.toString(),
          struct_asym_id: chain,
          start_residue_number: Number(elem.split('-')[0]),
          end_residue_number: Number(elem.split('-')[1]),
          color: { r: 23, g: 162, b: 184 },
        });
      });
      console.log(coloring);
      var theme = LiteMolPluginInstance.CustomTheme.createTheme(
        model.props.model,
        coloring
      );
      LiteMolPluginInstance.CustomTheme.applyTheme(
        plugin,
        'polymer-visual',
        theme
      );
    });
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
