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
