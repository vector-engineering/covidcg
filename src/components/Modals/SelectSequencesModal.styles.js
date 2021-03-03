import styled from 'styled-components';
import { lighten } from 'polished';

import Button from '../Buttons/Button';

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: flex-start;

  width: calc(100vw - 100px);
  height: calc(100vh - 100px);
`;

export const HeaderContainer = styled.div`
  position: fixed;
  top: 0;
  left: 0;
  width: 100%;
  flex-shrink: 0;
  background-color: #fff;

  height: 40px;

  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-start;

  border-bottom: 1px solid #ccc;
  box-shadow: 0px 0px 5px 1px #ddd;
`;

export const HeaderRow = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  height: 100%;
  width: 100%;
  padding: 3px 5px;
`;

export const HeaderButtons = styled.div`
  display: flex;
  flex-direction: row;
`;

export const TitleContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: flex-start;
  margin-left: 20px;

  .title {
    h2 {
      margin-bottom: 0px;
      margin-top: 0px;
    }
  }
  .spacer {
    flex-grow: 1;
  }
  .close-button {
  }
`;
export const Content = styled.div`
  height: 100%;

  font-size: 1em;
  font-weight: normal;

  display: flex;
  flex-direction: row;
  align-items: flex-start;
  flex-wrap: wrap;
`;

export const Column = styled.div`
  min-width: ${({ minWidth }) => minWidth}px;
  ${({ maxWidth }) => (maxWidth > 0 ? 'max-width: ' + maxWidth + 'px;' : '')}
  height: 100%;

  display: flex;
  flex-direction: column;
  align-items: stretch;
  overflow-y: scroll;
  ${({ grow }) => (grow ? 'flex-grow: 1;' : '')}
`;
Column.defaultProps = {
  minWidth: 300,
  grow: false,
  maxWidth: 0,
};

export const CancelButton = styled(Button)`
  background-color: #ddd;
  color: #000;
  background-image: none;
  font-size: 1rem;
  font-weight: normal;
  margin-right: 10px;

  &:hover,
  &:active {
    background-color: #eee;
  }
`;

export const ApplyButton = styled(Button)`
  background-image: none;
  font-size: 1rem;
  font-weight: normal;

  background-color: ${({ invalid }) => (invalid ? '#DDD' : '#28a745')};
  color: ${({ invalid }) => (invalid ? '#888' : '#FFF')};

  transition: 0.1s all ease-in-out;

  &:hover {
    background-color: ${({ invalid }) =>
      invalid ? '#DDD' : lighten(0.1, '#28a745')};
  }
`;
ApplyButton.defaultProps = {
  invalid: false,
};

export const InvalidText = styled.span`
  margin-reft: 15px;
  font-size: 1rem;
  font-weight: normal;
  line-height: normal;
  color: #dc3545;
`;

export const Overlay = styled.div`
  visibility: ${({ visible }) => (visible ? 'visible' : 'hidden')};
  position: absolute;
  top: 0;
  left: 0;
  width: 100%;
  height: 100%;
  z-index: 2;
  background-color: rgba(230, 230, 230, 0.8);

  display: flex;
  flex-direction: column;
  justify-content: center;
`;

export const ProgressContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: center;
`;

export const ProgressText = styled.span`
  margin-left: 20px;
  font: 14px 'Helvetica Neue', Helvetica, Arial, sans-serif;
  line-height: 1.4em;
  color: #4d4d4d;

  font-size: 2rem;
  font-weight: 500;
`;
