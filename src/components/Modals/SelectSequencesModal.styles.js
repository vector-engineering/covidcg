import styled from 'styled-components';
import { lighten } from 'polished';

import Button from '../Buttons/Button';

const FOLD_WIDTH = '1000px';

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

  display: grid;
  grid-template-rows: auto;
  grid-template-columns: repeat(3, 1fr);

  @media (max-width: ${FOLD_WIDTH}) {
    grid-template-columns: repeat(2, 1fr);
  }
`;

export const Column = styled.div`
  min-width: ${({ minWidth }) => minWidth}px;
  height: 100%;

  display: flex;
  flex-direction: column;
  align-items: stretch;
  overflow-y: scroll;

  @media (max-width: ${FOLD_WIDTH}) {
    grid-row: ${({ collapseRow }) => collapseRow};
    grid-column: ${({ collapseCol }) => collapseCol};
  }
`;
Column.defaultProps = {
  minWidth: 300,
  collapseRow: 1,
  collapseCol: 1,
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
