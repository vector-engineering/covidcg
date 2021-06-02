import LiteMol from 'litemol';

const Core = LiteMol.Core;
const Visualization = LiteMol.Visualization;
const Bootstrap = LiteMol.Bootstrap;
const Query = Core.Structure.Query;

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

function queryResiduesByIndices(indices) {
  return Query.Builder.build(() => compileResiduesByIndices(indices));
}

function compileResiduesByIndices(indices) {
  return (ctx) => {
    let { residues, chains } = ctx.structure.data,
      { authSeqNumber, atomStartIndex, atomEndIndex } = residues,
      { count, residueStartIndex, residueEndIndex } = chains,
      fragments = new Query.FragmentSeqBuilder(ctx);

    // console.log(
    //   seqNumber,
    //   authSeqNumber,
    //   insCode,
    //   atomStartIndex,
    //   atomEndIndex
    // );
    // console.log(entityId, count, residueStartIndex, residueEndIndex);

    const seqSource = authSeqNumber;

    for (let cI = 0; cI < count; cI++) {
      let i = residueStartIndex[cI],
        last = residueEndIndex[cI];

      let indicesI = 0;
      while (i < last && indicesI < indices.length - 1) {
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

export function createTheme(model, colorDef) {
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

export const colorHeatmap = ({ plugin, entries }) => {
  let model = plugin.selectEntities('model')[0];
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

  const theme = createTheme(model.props.model, coloring);

  // instead of "polymer-visual", "model" or any valid ref can be used: all "child" visuals will be colored.
  applyTheme(plugin, 'model', theme);
};

export const colorSequences = (plugin) => {
  let model = plugin.selectEntities('model')[0];
  if (!model) return;

  // hardcoded example
  const coloring = {
    base: { r: 255, g: 255, b: 255 },
    entries: [
      // {
      //   entity_id: '1',
      //   struct_asym_id: 'A',
      //   start_residue_number: 10,
      //   end_residue_number: 250,
      //   color: { r: 255, g: 128, b: 64 },
      // },
      // {
      //   entity_id: '1',
      //   struct_asym_id: 'A',
      //   start_residue_number: 300,
      //   end_residue_number: 600,
      //   color: { r: 64, g: 128, b: 255 },
      // },
      {
        indices: [...Array(100).keys()],
        color: { r: 64, g: 128, b: 255 },
      },
    ],
  };

  const theme = createTheme(model.props.model, coloring);

  // instead of "polymer-visual", "model" or any valid ref can be used: all "child" visuals will be colored.
  applyTheme(plugin, 'model', theme);
};
