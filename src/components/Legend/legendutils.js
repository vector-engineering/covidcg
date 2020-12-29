import {
  COORDINATE_MODES,
  DNA_OR_AA,
  appConfig,
  GROUP_SNV,
} from '../../constants/config';

export const sortLegendItems = (groupKey, dnaOrAa, coordinateMode, a, b) => {
  // If we're grouping by lineage or clade, then sort alphabetically
  // on the lineage/clade
  if (Object.keys(appConfig.group_cols).includes(groupKey)) {
    return a.group > b.group ? 1 : -1;
  } else if (groupKey === GROUP_SNV) {
    // If we're grouping by SNV, figure out whether we're in DNA or AA mode
    if (dnaOrAa === DNA_OR_AA.DNA) {
      // If we're grouping by DNA SNV, then sort by position
      return a.pos > b.pos ? 1 : -1;
    } else {
      // If we're grouping by AA SNV, determine if we're in gene/protein mode
      if (coordinateMode === COORDINATE_MODES.COORD_GENE) {
        // If the same gene, then sort by position
        if (a.gene === b.gene) {
          return a.pos > b.pos ? 1 : -1;
        }
        // Otherwise, sort by gene
        else {
          return a.gene > b.gene ? 1 : -1;
        }
      } else if (coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
        // If the same protein, then sort by position
        if (a.protein === b.protein) {
          return a.pos > b.pos ? 1 : -1;
        }
        // Otherwise, sort by protein
        else {
          return a.protein.toLowerCase() > b.protein.toLowerCase() ? 1 : -1;
        }
      }
    }
  }
};
