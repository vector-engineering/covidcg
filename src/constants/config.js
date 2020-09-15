const LOCAL_COUNTS = 'local';
const GLOBAL_COUNTS = 'global';
const GROUP_COUNTS = 'group';
export const LOW_FREQ_FILTER_TYPES = {
  LOCAL_COUNTS,
  GLOBAL_COUNTS,
  GROUP_COUNTS,
};

const GROUP_LINEAGE = 'lineage';
const GROUP_CLADE = 'clade';
const GROUP_SNV = 'snv';
export const GROUP_KEYS = {
  GROUP_LINEAGE,
  GROUP_CLADE,
  GROUP_SNV,
};

const DNA = 'DNA';
const AA = 'AA';
export const DNA_OR_AA = {
  DNA,
  AA,
};

const COORD_GENE = 'gene';
const COORD_PROTEIN = 'protein';
const COORD_PRIMER = 'primer';
const COORD_CUSTOM = 'custom';
const COORD_SEQUENCE = 'sequence';
export const COORDINATE_MODES = {
  COORD_GENE,
  COORD_PROTEIN,
  COORD_PRIMER,
  COORD_CUSTOM,
  COORD_SEQUENCE,
};
