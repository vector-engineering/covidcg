import LiteMol from 'litemol';

import { LITEMOL_STYLES } from '../../constants/defs';

const Core = LiteMol.Core;
const Visualization = LiteMol.Visualization;
const Bootstrap = LiteMol.Bootstrap;
const Query = Core.Structure.Query;
const Tree = Bootstrap.Tree;
const Transformer = Bootstrap.Entity.Transformer;

// const Builder = Query.Builder;
// const Compiler = Query.Compiler;
// class OptimizedId {
//   columns;
//   isSatisfied(i) {
//     for (let c of this.columns) {
//       if (c.value !== c.array[i]) return false;
//     }
//     return true;
//   }
//   constructor(id, arrays) {
//     this.columns = [];
//     for (let key of Object.keys(id)) {
//       if (id[key] !== void 0 && !!arrays[key]) {
//         this.columns.push({ value: id[key], array: arrays[key] });
//       }
//     }
//   }
// }

function queryAtomsByChains(chainIds) {
  return Query.Builder.build(() => compileAtomsByChains(chainIds));
}

function compileAtomsByChains(chainIds) {
  return (ctx) => {
    let { chains } = ctx.structure.data,
      { authAsymId, count, atomStartIndex, atomEndIndex } = chains,
      fragments = new Query.FragmentSeqBuilder(ctx);

    for (let chainI = 0; chainI < count; chainI++) {
      const chainAsymId = authAsymId[chainI];

      if (chainIds.includes(chainAsymId)) {
        fragments.add(
          Query.Fragment.ofIndexRange(
            ctx,
            atomStartIndex[chainI],
            atomEndIndex[chainI]
          )
        );
      }
    }

    return fragments.getSeq();
  };
}

function queryResiduesByIndices(indices) {
  return Query.Builder.build(() => compileResiduesByIndices(indices));
}

function compileResiduesByIndices(indices) {
  return (ctx) => {
    let { residues, chains } = ctx.structure.data,
      { authSeqNumber, atomStartIndex, atomEndIndex } = residues,
      { count, residueStartIndex, residueEndIndex } = chains,
      fragments = new Query.FragmentSeqBuilder(ctx);

    // console.log(authSeqNumber, atomStartIndex, atomEndIndex);
    // console.log(count, residueStartIndex, residueEndIndex);

    const seqSource = authSeqNumber;

    // For each chain
    for (let chainI = 0; chainI < count; chainI++) {
      // Get the residue indices for this chain
      let i = residueStartIndex[chainI],
        lastResidueIndex = residueEndIndex[chainI];

      let indicesI = 0;
      while (i < lastResidueIndex && indicesI <= indices.length - 1) {
        // If the current residue index is below the
        // active query residue index, then increment the current
        // residue counter
        if (seqSource[i] < indices[indicesI]) {
          i += 1;
        }
        // If we have a match, then add fragment to query
        // and increment both our current residue counter
        // and the query residue index
        else if (seqSource[i] == indices[indicesI]) {
          if (ctx.hasRange(atomStartIndex[i], atomEndIndex[i])) {
            fragments.add(
              Query.Fragment.ofIndexRange(
                ctx,
                atomStartIndex[i],
                atomEndIndex[i]
              )
            );
          }
          i += 1;
          indicesI += 1;
        }
        // If the current residue index is above
        // the active query residue index, then we skipped it
        // (most likely the query residue index doesn't exist)
        // Increment the query residue index and try again
        else if (seqSource[i] > indices[indicesI]) {
          indicesI += 1;
        }
      }
    }

    return fragments.getSeq();
  };
}

// function sequence(entityId, asymId, startId, endId) {
//   return Query.Builder.build(() =>
//     compileSequence(entityId, asymId, startId, endId)
//   );
// }

// function compileSequence(seqEntityId, seqAsymId, start, end) {
//   return (ctx) => {
//     let { residues, chains } = ctx.structure.data,
//       { seqNumber, authSeqNumber, insCode, atomStartIndex, atomEndIndex } =
//         residues,
//       { entityId, count, residueStartIndex, residueEndIndex } = chains,
//       fragments = new Query.FragmentSeqBuilder(ctx);

//     let parent = ctx.structure.parent,
//       { sourceChainIndex } = chains,
//       isComputed = parent && sourceChainIndex;

//     let targetAsymId =
//       typeof seqAsymId === 'string' ? { asymId: seqAsymId } : seqAsymId;
//     let optTargetAsymId = new OptimizedId(
//       targetAsymId,
//       isComputed ? parent.data.chains : chains
//     );

//     console.log(
//       seqNumber,
//       authSeqNumber,
//       insCode,
//       atomStartIndex,
//       atomEndIndex
//     );
//     console.log(entityId, count, residueStartIndex, residueEndIndex);

//     // for (let i of indices) {
//     //     if (!ctx.hasRange(atomStartIndex[i], atomEndIndex[i])) continue;
//     //     fragments.add(Fragment.ofIndexRange(ctx, atomStartIndex[i], atomEndIndex[i]));
//     // }

//     const isAuth = typeof targetAsymId.authAsymId === 'string';

//     const seqSource = isAuth ? authSeqNumber : seqNumber;
//     const startSeqNumber = isAuth ? start.authSeqNumber : start.seqNumber;
//     const endSeqNumber = isAuth ? end.authSeqNumber : end.seqNumber;

//     //optAsymId.isSatisfied();

//     for (let cI = 0; cI < count; cI++) {
//       // if (
//       //   (!!seqEntityId && entityId[cI] !== seqEntityId) ||
//       //   !optTargetAsymId.isSatisfied(isComputed ? sourceChainIndex[cI] : cI)
//       // ) {
//       //   continue;
//       // }

//       let i = residueStartIndex[cI],
//         last = residueEndIndex[cI],
//         startIndex = -1,
//         endIndex = -1;
//       for (; i < last; i++) {
//         if (seqSource[i] >= startSeqNumber && seqSource[i] <= endSeqNumber) {
//           if (!!start.insCode && insCode[i] !== start.insCode) continue;
//           startIndex = i;
//           break;
//         }
//       }

//       if (i < 0 || i === last) continue;

//       for (i = startIndex; i < last; i++) {
//         if (seqSource[i] >= endSeqNumber) {
//           if (
//             !!end.insCode &&
//             seqSource[i] === endSeqNumber &&
//             insCode[i] !== end.insCode
//           )
//             continue;
//           break;
//         }
//       }

//       endIndex = i;

//       if (ctx.hasRange(atomStartIndex[startIndex], atomEndIndex[endIndex])) {
//         fragments.add(
//           Query.Fragment.ofIndexRange(
//             ctx,
//             atomStartIndex[startIndex],
//             atomEndIndex[endIndex]
//           )
//         );
//       }
//     }

//     return fragments.getSeq();
//   };
// }

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

export function createTheme(model, colorDef, ignoreChains) {
  const mapper = new ColorMapper();
  mapper.addColor(colorDef.base);
  const map = new Uint8Array(model.data.atoms.count);

  for (const e of colorDef.entries) {
    // const query = Query.sequence(
    //   e.entity_id.toString(),
    //   e.struct_asym_id,
    //   { seqNumber: e.start_residue_number },
    //   { seqNumber: e.end_residue_number }
    // ).compile();
    // const query = sequence(
    //   e.entity_id.toString(),
    //   e.struct_asym_id,
    //   { seqNumber: e.start_residue_number },
    //   { seqNumber: e.end_residue_number }
    // ).compile();
    const query = queryResiduesByIndices(e.indices).compile();
    const colorIndex = mapper.addColor(e.color);
    for (const f of query(model.queryContext).fragments) {
      for (const a of f.atomIndices) {
        map[a] = colorIndex;
      }
    }
  }

  const ignoreChainColor = { r: 150, g: 150, b: 150 };
  const colorIndex = mapper.addColor(ignoreChainColor);
  // console.log(ignoreChains);
  if (ignoreChains.length > 0) {
    const query = queryAtomsByChains(ignoreChains).compile();
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
      .filter((node) => node.parent.props.label == 'Polymer')
  );

  for (const v of visuals) {
    plugin.command(Bootstrap.Command.Visual.UpdateBasicTheme, {
      visual: v,
      theme,
    });
  }
}

export const colorHeatmap = ({ plugin, entries, ref, ignoreChains }) => {
  const model = plugin.selectEntities(ref)[0];
  if (!model) return;

  // hardcoded example
  const coloring = {
    base: { r: 255, g: 255, b: 255 },
    // [{
    //   indices: [1, 2, 3, 4, 5],
    //   color: { r: 64, g: 128, b: 255 },
    // }, ... ]
    entries,
  };

  const theme = createTheme(model.props.model, coloring, ignoreChains);

  // instead of "polymer-visual", "model" or any valid ref can be used: all "child" visuals will be colored.
  applyTheme(plugin, ref, theme);
};

// This might help... https://github.com/dsehnal/LiteMol/blob/2ce0190a9b369841c1c3ee7322c2f5dda3e7800e/examples/CustomControls/src/Representation.ts
export const getMoleculeAssemblies = ({ plugin }) => {
  const molecule = plugin.selectEntities('molecule')[0];
  if (!molecule) return [];

  const assemblies =
    molecule.props.molecule.models[0].data.assemblyInfo.assemblies;
  return assemblies.map((asm) => asm.name);
};

export const getMoleculeEntities = ({ plugin }) => {
  const molecule = plugin.selectEntities('molecule')[0];
  if (!molecule) return [];

  // See: https://github.com/dsehnal/LiteMol/blob/2ce0190a9b369841c1c3ee7322c2f5dda3e7800e/src/lib/Core/lib/CIFTools.js
  //      https://github.com/dsehnal/LiteMol/blob/2ce0190a9b369841c1c3ee7322c2f5dda3e7800e/src/lib/Core/lib/CIFTools.d.ts
  const categoryMap =
    molecule.parent.props.dictionary.dataBlocks[0].categoryMap;

  const entity = categoryMap.get('_entity').toJSON().rows;
  const entityPoly = categoryMap.get('_entity_poly').toJSON().rows;

  // console.log(entity, entityPoly);

  const entityObjs = entity
    .map((e) => {
      e = Object.assign({}, e);

      // Skip over non-polymers. We won't be coloring het groups anyways
      if (e.type !== 'polymer') {
        return null;
      }

      // Find the corresponding chains and polymer type from the entityPoly rows
      const poly = entityPoly.find((p) => p.entity_id === e.id);

      // If no poly object found, break out
      if (poly === undefined) {
        return null;
      }

      // Attach chain and polymer type information
      // and also the sequence I guess
      e.chains = poly.pdbx_strand_id.split(',');
      e.poly_type = poly.type;
      // Check for key ownership just in case... maybe this polymer is DNA or something
      if (
        Object.prototype.hasOwnProperty.call(poly, 'pdbx_seq_one_letter_code')
      ) {
        e.seq = poly.pdbx_seq_one_letter_code;
      }

      // Apply heatmap by default?
      e.checked = true;

      // Disable some entities by default
      const descriptionBlacklist = [
        // Antibody-like
        /heavy\schain/gi,
        /light\schain/gi,
        /[FNM]ab/g,
        // polynucleicacids
        /[DR]NA$/g,
      ];
      e.checked = descriptionBlacklist.every((regex) => {
        return e.pdbx_description.search(regex) === -1;
      });
      // Disable ribonucleotide types
      const typeBlacklist = ['polyribonucleotide'];
      e.checked =
        e.checked &&
        typeBlacklist.every((t) => {
          return e.poly_type !== t;
        });

      return e;
    })
    .filter((e) => e !== null);

  // console.log(entityObjs);

  return entityObjs;
};

const selectionColors = Bootstrap.Immutable.Map()
  .set('Uniform', Visualization.Color.fromHex(0xaaaaaa))
  .set('Selection', Visualization.Theme.Default.SelectionColor)
  .set('Highlight', Visualization.Theme.Default.HighlightColor);

const polymerSurfaceStyle = {
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

const polymerCartoonStyle = {
  type: 'Cartoons',
  params: {
    detail: 'Very High',
    showDirectionCone: false,
  },
  theme: {
    template: Bootstrap.Visualization.Molecule.Default.UniformThemeTemplate,
    colors: selectionColors,
    transparency: { alpha: 1.0 },
  },
};

// const ligandStyle = {
//   type: 'BallsAndSticks',
//     params: {
//       useVDW: true,
//       vdwScaling: 0.25,
//       bondRadius: 0.13,
//       detail: 'Automatic'
//     },
//     theme: {
//       template: Bootstrap.Visualization.Molecule.Default.ElementSymbolThemeTemplate,
//       colors: Bootstrap.Visualization.Molecule.Default.ElementSymbolThemeTemplate.colors,
//       transparency: { alpha: 1.0 }
//     },
// };

// const waterStyle = {
//   type: 'BallsAndSticks',
//   params: {
//     useVDW: false,
//     atomRadius: 0.23,
//     bondRadius: 0.09,
//     detail: 'Automatic'
//   },
//   theme: {
//     template: Bootstrap.Visualization.Molecule.Default.ElementSymbolThemeTemplate,
//     colors: Bootstrap.Visualization.Molecule.Default.ElementSymbolThemeTemplate.colors,
//     transparency: { alpha: 0.25 }
//   },
// };

export const CreateMacromoleculeVisual = Tree.Transformer.action(
  {
    id: 'molecule-create-macromolecule-visual',
    name: 'Macromolecule Visual',
    description:
      'Create a visual of a molecule that is split into polymer, HET, and water parts.',
    from: [
      Bootstrap.Entity.Molecule.Selection,
      Bootstrap.Entity.Molecule.Model,
    ],
    to: [Bootstrap.Entity.Action],
    validateParams: (p) =>
      !p.polymer && !p.het && !p.water
        ? ['Select at least one component']
        : void 0,
    // eslint-disable-next-line no-unused-vars
    defaultParams: (_ctx) => ({
      polymer: true,
      het: true,
      water: true,
      style: LITEMOL_STYLES.SURFACE,
    }),
  },
  (_context, a, t) => {
    let g = Tree.Transform.build().add(
      a,
      Bootstrap.Entity.Transformer.Basic.CreateGroup,
      { label: 'Group', description: 'Macromolecule' },
      { ref: t.params.groupRef }
    );

    if (t.params.polymer) {
      let polymerStyle;
      if (t.params.style === LITEMOL_STYLES.SURFACE) {
        polymerStyle = polymerSurfaceStyle;
      } else if (t.params.style === LITEMOL_STYLES.CARTOON) {
        polymerStyle = polymerCartoonStyle;
      }

      g.then(
        Transformer.Molecule.CreateSelectionFromQuery,
        {
          query: Core.Structure.Query.nonHetPolymer(),
          name: 'Polymer',
          silent: true,
        },
        { isBinding: true }
      ).then(
        Transformer.Molecule.CreateVisual,
        {
          style: polymerStyle,
        },
        { ref: t.params.polymerRef }
      );
    }

    if (t.params.het) {
      g.then(
        Transformer.Molecule.CreateSelectionFromQuery,
        { query: Core.Structure.Query.hetGroups(), name: 'HET', silent: true },
        { isBinding: true }
      ).then(
        Transformer.Molecule.CreateVisual,
        {
          style:
            Bootstrap.Visualization.Molecule.Default.ForType.get(
              'BallsAndSticks'
            ),
        },
        { ref: t.params.hetRef }
      );
    }

    if (t.params.water) {
      let style = {
        type: 'BallsAndSticks',
        params: {
          useVDW: false,
          atomRadius: 0.23,
          bondRadius: 0.09,
          detail: 'Automatic',
        },
        theme: {
          template: Visualization.Molecule.Default.ElementSymbolThemeTemplate,
          colors:
            Bootstrap.Visualization.Molecule.Default.ElementSymbolThemeTemplate
              .colors,
          transparency: { alpha: 0.25 },
        },
      };

      g.then(
        Transformer.Molecule.CreateSelectionFromQuery,
        {
          query: Core.Structure.Query.entities({ type: 'water' }),
          name: 'Water',
          silent: true,
        },
        { isBinding: true }
      ).then(
        Transformer.Molecule.CreateVisual,
        { style },
        { ref: t.params.waterRef }
      );
    }
    return g;
  }
);

export const LoadLitemolModel = ({
  plugin,
  pdbId,
  proteinStyle,
  useAssembly,
  onLoad,
}) => {
  plugin.clear();
  pdbId = pdbId.toLowerCase();

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

    vizAction = vizAction.then(CreateMacromoleculeVisual, {
      polymer: true,
      het: true,
      water: false,
      style: proteinStyle,
    });

    onLoad({
      vizAction,
      assemblies,
      entities,
      activeAssembly: useAssembly,
    });
  });

  return modelAction;
};
