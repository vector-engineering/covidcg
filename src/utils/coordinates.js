import { config } from '../config';

export function coordsToText(coordinateRanges, showSegment) {
  if (showSegment === undefined) {
    showSegment = config.segments.length > 1;
  }

  return coordinateRanges
    .map((range) => {
      if (showSegment) {
        return `${range[0]}:${range[1]}..${range[2]}`;
      } else {
        return `${range[1]}..${range[2]}`;
      }
    })
    .join(';');
}

export function textToCoords(text) {
  return text.split(';').map((range) => {
    // If the range doesn't contain a segment, then default to first segment
    let segment;
    if (range.split(':').length === 1) {
      segment = config.segments[0]; // Default to first segment
    } else {
      segment = range.split(':')[0];
      range = range.split(':')[1];
    }

    // Extract start/end
    //const [start, end] = range.split('..').map((coord) => parseInt(coord));
    const [start, end] = range.split('..');

    return [segment, start, end];
  });
}

// Residue coordinates more straightforward -- no segment information
export function residueCoordsToText(residueCoordinates) {
  return residueCoordinates
    .map((coord) => {
      return `${coord[0]}..${coord[1]}`;
    })
    .join(';');
}

export function textToResidueCoords(text) {
  return text.split(';').map((range) => {
    const [start, end] = range.split('..').map((coord) => parseInt(coord));
    return [start, end];
  });
}
