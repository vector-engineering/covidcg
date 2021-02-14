// The configuration YAML file (defined by the environment variable CONFIGFILE
// when the app build is launched) is injected in place of the string "CG_CONFIG"
// as a JSON object, via webpack.DefinePlugin
const config = CG_CONFIG; // eslint-disable-line no-undef

export const appConfig = Object.assign({}, config);

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

export const MIN_DATE = 1576368000000;
