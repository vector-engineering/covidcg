import styled from 'styled-components';
import { lighten, transparentize } from 'polished';

export const LegendList = styled.div`
  display: flex;
  flex-direction: row;
  flex-wrap: wrap;
  border-radius: 2px;
  padding: 6px 3px;
`;

export const LegendItem = styled.div`
  display: flex;
  align-items: center;
  padding: 0px 6px;
  user-select: none;

  border: ${({ hovered, selected }) => {
    if (hovered) {
      return '1px solid #666';
    } else if (selected) {
      return '1px solid #000';
    } else {
      return 'none';
    }
  }};
  border-radius: 3px;
  margin: ${({ hovered, selected }) => {
    if (hovered || selected) {
      return '0px 4px 3px 1px';
    } else {
      return '1px 5px 4px 2px';
    }
  }};
  background-color: ${({ color, hovered, selected }) => {
    if (hovered) {
      return lighten(0.1, color);
    } else if (selected === null) {
      return color;
    } else if (selected) {
      return color;
    } else {
      return transparentize(0.7, color);
    }
  }};

  font-size: 12px;
  font-weight: 500;
  color: ${({ textColor, hovered, selected }) => {
    if (hovered) {
      return transparentize(0.1, textColor);
    } else if (selected === null) {
      return textColor;
    } else if (selected) {
      return textColor;
    } else {
      // Make the blacks a lot lighter
      return textColor === '#fff'
        ? transparentize(0.5, textColor)
        : transparentize(0.8, textColor);
    }
  }};
`;

LegendItem.defaultProps = {
  hovered: false,
  selected: null,
  color: '#ccc',
  textColor: '#fff',
};

export const ColorCircle = styled.div`
  background-color: ${({ color, selected }) => {
    if (selected === null || selected) {
      return color;
    } else {
      return '#CCC';
    }
  }};
  width: 12px;
  height: 12px;
  border-radius: 50%;
  margin-right: 4px;
`;

ColorCircle.defaultProps = {
  color: '#000',
  selected: null,
};
