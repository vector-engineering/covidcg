import config from '../../config.yaml';

export const appConfig = Object.assign({}, config);
export const GROUP_COLS = config.group_cols[config.ingest_strategy];
export const METADATA_COLS = config.metadata_cols[config.ingest_strategy];

const LOCAL_COUNTS = 'local';
const GLOBAL_COUNTS = 'global';
const GROUP_COUNTS = 'group';
export const LOW_FREQ_FILTER_TYPES = {
  LOCAL_COUNTS,
  GLOBAL_COUNTS,
  GROUP_COUNTS,
};

export const GROUP_SNV = 'snv';

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
