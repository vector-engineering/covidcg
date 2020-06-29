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
