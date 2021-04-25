import LiteMol from 'litemol';

const Core = LiteMol.Core;
const Visualization = LiteMol.Visualization;
const Bootstrap = LiteMol.Bootstrap;
const Q = Core.Structure.Query;

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
    if (this.map.has(id)) return this.map.get(id);
    const index = this.uniqueColors.length;
    this.uniqueColors.push(
      Visualization.Color.fromRgb(color.r, color.g, color.b)
    );
    this.map.set(id, index);
    return index;
  }
}

export function createTheme(model, colorDef) {
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
  return Visualization.Theme.createMapping(mapping, { colors, isSticky: true });
}

export function applyTheme(plugin, modelRef, theme) {
  const visuals = plugin.selectEntities(
    Bootstrap.Tree.Selection.byRef(modelRef)
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

export const colorSequences = (plugin) => {
  let model = plugin.selectEntities('model')[0];
  if (!model) return;

  // hardcoded example
  const coloring = {
    base: { r: 255, g: 255, b: 255 },
    entries: [
      {
        entity_id: '1',
        struct_asym_id: 'A',
        start_residue_number: 10,
        end_residue_number: 250,
        color: { r: 255, g: 128, b: 64 },
      },
      {
        entity_id: '1',
        struct_asym_id: 'A',
        start_residue_number: 300,
        end_residue_number: 600,
        color: { r: 64, g: 128, b: 255 },
      },
    ],
  };

  const theme = createTheme(model.props.model, coloring);

  // instead of "polymer-visual", "model" or any valid ref can be used: all "child" visuals will be colored.
  applyTheme(plugin, 'model', theme);
};
