import styled from 'styled-components';

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
