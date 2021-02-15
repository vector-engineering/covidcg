import {
  COORDINATE_MODES,
  DNA_OR_AA,
  GROUP_SNV,
} from '../../constants/defs.json';
import { config } from '../../config';

export const sortLegendItems = (groupKey, dnaOrAa, coordinateMode, a, b) => {
  // If we're grouping by lineage or clade, then sort alphabetically
  // on the lineage/clade
  if (Object.keys(config.group_cols).includes(groupKey)) {
    return a.group > b.group ? 1 : -1;
  } else if (groupKey === GROUP_SNV) {
    // If we're grouping by SNV, figure out whether we're in DNA or AA mode
    if (dnaOrAa === DNA_OR_AA.DNA) {
      // If we're grouping by DNA SNV, then sort by position
      return a.pos > b.pos ? 1 : -1;
    } else {
      // If we're grouping by AA SNV, determine if we're in gene/protein mode
      if (
        coordinateMode === COORDINATE_MODES.COORD_GENE ||
        coordinateMode === COORDINATE_MODES.COORD_PROTEIN
      ) {
        // If the same gene/protein, then sort by position
        if (a.name === b.name) {
          return a.pos > b.pos ? 1 : -1;
        }
        // Otherwise, sort by name of the gene/protein
        else {
          return a.name.toLowerCase() > b.name.toLowerCase() ? 1 : -1;
        }
      }
    }
  }
};
