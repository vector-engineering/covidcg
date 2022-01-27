import { COMPARE_COLORS } from '../constants/defs.json';

export const warmColors = [
  // autumn palette, [0, 1, 0.15]
  // '#ff0000',
  // '#ff2600',
  // '#ff4c00',
  // '#ff7300',
  // '#ff9900',
  // '#ffc000',
  // '#ffe600',
  // SHUFFLED:
  '#ff0000',
  '#ff7300',
  '#ff2600',
  '#ff9900',
  '#ff4c00',
  '#ffc000',
  '#ffe600',
];

export const coolColors = [
  // bwr palette, [0, 0.4, 0.07]
  // '#0000ff',
  // '#2222ff',
  // '#4646ff',
  // '#6969ff',
  // '#8e8eff',
  // '#b2b2ff',
  // cool palette, [0.2, 0.8, 0.05]
  // '#33ccff',
  // '#40bfff',
  // '#4cb3ff',
  // '#59a6ff',
  // '#6699ff',
  // '#738cff',
  // '#7f80ff',
  // '#8c72ff',
  // '#9966ff',
  // '#a659ff',
  // '#b34cff',
  // '#bf40ff',
  // '#cc32ff',
  // SHUFFLED:
  '#0000ff',
  '#6969ff',
  '#2222ff',
  '#8e8eff',
  '#4646ff',
  '#b2b2ff',
  '#33ccff',
  '#7f80ff',
  '#40bfff',
  '#8c72ff',
  '#4cb3ff',
  '#9966ff',
  '#59a6ff',
  '#a659ff',
  '#6699ff',
  '#b34cff',
  '#738cff',
  '#bf40ff',
];

// from https://personal.sron.nl/~pault/#sec:qualitative
// 'vibrant', then 'muted', then 'light'
export const mutationColorArray = [
  '#0077bb',
  '#33bbee',
  '#009988',
  '#ee7733',
  '#cc3311',
  '#ee3377',

  '#332288',
  '#88ccee',
  '#44aa99',
  '#117733',
  '#999933',
  '#ddcc77',
  '#cc6677',
  '#882255',
  '#aa4499',

  '#77AADD',
  '#99DDFF',
  '#44BB99',
  '#BBCC33',
  '#AAAA00',
  '#EEDD88',
  '#EE8866',
  '#FFAABB',
  '#DDDDDD',
];

// from https://personal.sron.nl/~pault/#sec:qualitative
// 'vibrant' scheme
export const cladeColorArray = [
  '#0077bb',
  '#33bbee',
  '#009988',
  '#ee7733',
  '#cc3311',
  '#ee3377',
];

export const snapGeneNTColors = {
  A: '#6afb60',
  C: '#9acaff',
  G: '#fee400',
  T: '#f096ac',
};

export const snapGeneHighlightColors = {};
snapGeneHighlightColors[COMPARE_COLORS.COMPARE_COLOR_YELLOW] = '#ffe928';
snapGeneHighlightColors[COMPARE_COLORS.COMPARE_COLOR_GREEN] = '#beed64';
snapGeneHighlightColors[COMPARE_COLORS.COMPARE_COLOR_BLUE] = '#a9dfff';
snapGeneHighlightColors[COMPARE_COLORS.COMPARE_COLOR_PINK] = '#fdc2c3';
snapGeneHighlightColors[COMPARE_COLORS.COMPARE_COLOR_PURPLE] = '#d5ccff';
snapGeneHighlightColors[COMPARE_COLORS.COMPARE_COLOR_ORANGE] = '#ffd900'; // more of a yellow than an orange...
snapGeneHighlightColors[COMPARE_COLORS.COMPARE_COLOR_GRAY] = '#d1d1d1';

export const shingAAColors = {
  R: '#E60606',
  K: '#C64200',
  Q: '#FF6600',
  N: '#CCFF99',
  E: '#FFCC00',
  D: '#FFCC99',
  H: '#FFFF99',
  P: '#FFFF00',
  Y: '#CCFFCC',
  W: '#00FF00',
  S: '#FF9900',
  T: '#00FF99',
  G: '#CC99FF',
  A: '#CCFFFF',
  M: '#99CCFF',
  C: '#00FFFF',
  F: '#00CCFF',
  L: '#3366FF',
  V: '#0000FF',
  I: '#000080',
  _: '#FF0000',
};

// ClustalX AA Colors
// (by properties + conservation)
// http://www.jalview.org/help/html/colourSchemes/clustal.html
export const clustalXAAColors = {
  // Hydrophobic (Blue)
  A: '#809df0',
  I: '#809df0',
  L: '#809df0',
  M: '#809df0',
  F: '#809df0',
  W: '#809df0',
  V: '#809df0',
  // Positive charge (Red)
  K: '#ed000a',
  R: '#ed000a',
  // Negative charge (Magenta)
  D: '#be38bf',
  E: '#be38bf',
  // Polar (Green)
  N: '#29c417',
  Q: '#29c417',
  S: '#29c417',
  T: '#29c417',
  // Cysteins (Pink)
  C: '#ee7d80',
  // Glycines (Orange)
  G: '#ef8f48',
  // Prolines (Yellow)
  P: '#c1c204',
  // Aromatics (Cyan)
  H: '#23a6a4',
  Y: '#23a6a4',
  // STOP
  _: '#FF0000',
};

// Zappo Color Scheme (physico-chemical properties)
// From SnapGene
export const zappoAAColors = {
  A: '#fc8184',
  C: '#ffd900',
  D: '#eb3840',
  E: '#eb3840',
  F: '#fc7400',
  G: '#b858be',
  H: '#799af1',
  I: '#fc8184',
  K: '#799af1',
  L: '#fc8184',
  M: '#fc8184',
  N: '#25be00',
  P: '#b858be',
  Q: '#25be00',
  R: '#799af1',
  S: '#25be00',
  T: '#25be00',
  V: '#fc8184',
  W: '#fc7400',
  Y: '#fc7400',
  // STOP
  _: '#FF0000',
};

// Zhao and London (transmembrane-tendency)
// From SnapGene
export const transmembraneAAColors = {
  A: '#c85672',
  C: '#ba6288',
  D: '#7c96eb',
  E: '#8490df',
  F: '#eb3840',
  G: '#bc6085',
  H: '#a276ae',
  I: '#eb3840',
  K: '#799af1',
  L: '#e83b45',
  M: '#de4352',
  N: '#9e79b3',
  P: '#a276ae',
  Q: '#997dbb',
  R: '#8a8ad4',
  S: '#b56690',
  T: '#ba6289',
  V: '#e83b45',
  W: '#e1404e',
  Y: '#ca546f',
  // STOP
  _: '#FF0000',
};

export const reds = [
  '#FFF5F0',
  '#FEF1EB',
  '#FEEEE6',
  '#FEEAE1',
  '#FEE7DC',
  '#FEE3D7',
  '#FEE0D2',
  '#FDDACB',
  '#FDD4C3',
  '#FDCEBB',
  '#FCC8B3',
  '#FCC2AB',
  '#FCBCA3',
  '#FCB59B',
  '#FCAF93',
  '#FCA88B',
  '#FCA184',
  '#FC9B7C',
  '#FC9474',
  '#FB8D6D',
  '#FB8767',
];
