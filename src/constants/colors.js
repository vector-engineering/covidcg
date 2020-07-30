export const incrementColor = function (color, step) {
  var colorToInt = parseInt(color.substr(1), 16), // Convert HEX color to integer
    nstep = parseInt(step); // Convert step to integer
  if (!isNaN(colorToInt) && !isNaN(nstep)) {
    // Make sure that color has been converted to integer
    colorToInt += nstep; // Increment integer with step
    var ncolor = colorToInt.toString(16); // Convert back integer to HEX
    ncolor = '#' + new Array(7 - ncolor.length).join(0) + ncolor; // Left pad "0" to make HEX look like a color
    if (/^#[0-9a-f]{6}$/i.test(ncolor)) {
      // Make sure that HEX is a valid color
      return ncolor;
    }
  }
  return color;
};

export const warmColors = {
  //blood red
  colors: [
    '#6d0019',
    '#821b2a',
    '#96303c',
    '#ab434f',
    '#c05662',
    '#d56a77',
    '#ea7d8c',
    '#ff91a1',
  ],
  child: {
    // soft red
    colors: [
      '#ff5e5e',
      '#f45455',
      '#e9494c',
      '#de3e42',
      '#d4333a',
      '#c92731',
      '#be1828',
      '#b30020',
    ],
    child: {
      //orange
      colors: [
        '#ff8214',
        '#ff8f21',
        '#fe9c2d',
        '#fea83a',
        '#feb348',
        '#febe56',
        '#fec964',
        '#ffd373',
      ],
      // pink
      child: {
        colors: [
          '#ff1290',
          '#f40f88',
          '#ea0c81',
          '#df0979',
          '#d50672',
          '#cb046b',
          '#c00264',
          '#b6005d',
        ],
      },
    },
  },
};

export const coolColors = {
  // dark blue
  colors: [
    '#004c6d',
    '#006386',
    '#007b9f',
    '#0094b6',
    '#00aecb',
    '#00c9df',
    '#00e4f0',
    '#00ffff',
  ],
  child: {
    // light blue
    colors: [
      '#00666d',
      '#007a81',
      '#008f96',
      '#00a4ab',
      '#00bac0',
      '#00d1d5',
      '#00e8ea',
      '#00ffff',
    ],
    child: {
      //blue green
      colors: [
        '#006d50',
        '#038060',
        '#069471',
        '#08a982',
        '#0abe94',
        '#0ad3a7',
        '#07e9ba',
        '#00ffce',
      ],
      child: {
        // green
        colors: [
          '#096d0c',
          '#1b8018',
          '#2a9424',
          '#39a930',
          '#47be3c',
          '#55d348',
          '#63e954',
          '#71ff60',
        ],
      },
    },
  },
};

export const snpColorArray = [
  '#fb6b1d',
  '#e83b3b',
  '#831c5d',
  '#c32454',
  '#f04f78',
  '#f68181',
  '#fca790',
  '#e3c896',
  '#ab947a',
  '#966c6c',
  '#625565',
  '#0b5e65',
  '#0b8a8f',
  '#1ebc73',
  '#91db69',
  '#fbff86',
  '#fbb954',
  '#cd683d',
  '#9e4539',
  '#7a3045',
  '#6b3e75',
  '#905ea9',
  '#a884f3',
  '#eaaded',
  '#8fd3ff',
  '#4d9be6',
  '#4d65b4',
  '#484a77',
  '#30e1b9',
  '#8ff8e2',
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
  '#bbbbbb',
];

export const snapGeneNTColors = {
  A: '#6afb60',
  C: '#9acaff',
  G: '#fee400',
  T: '#f096ac',
};

export const snapGeneHighlightColors = {
  yellow: '#ffe928',
  green: '#beed64',
  blue: '#a9dfff',
  pink: '#fdc2c3',
  purple: '#d5ccff',
  orange: '#ffd900', // more of a yellow than an orange...
  gray: '#d1d1d1',
};

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
