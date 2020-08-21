import {
  warmColors,
  coolColors,
  snpColorArray,
  cladeColorArray,
} from '../constants/colors';
import _ from 'underscore';

let coolColorInd = 0;
let warmColorInd = 0;
export const getLineageColor = _.memoize((group) => {
  let color;
  if (group.charAt(0) === 'A') {
    color = warmColors[warmColorInd++];

    if (warmColorInd === warmColors.length) {
      warmColorInd = 0;
    }
  } else if (group.charAt(0) === 'B') {
    color = coolColors[coolColorInd++];

    if (coolColorInd === coolColors.length) {
      coolColorInd = 0;
    }
  }
  return color;
});

let snvColorInd = 0;
export const getSnvColor = _.memoize(() => {
  const color = snpColorArray[snvColorInd++];

  if (snvColorInd === snpColorArray.length) {
    snvColorInd = 0;
  }

  return color;
});

let cladeColorInd = 0;
export const getCladeColor = _.memoize(() => {
  const color = cladeColorArray[cladeColorInd++];
  // If we're at the end, then loop back to the beginning
  if (cladeColorInd === cladeColorArray.length) {
    cladeColorInd = 0;
  }

  return color;
});
