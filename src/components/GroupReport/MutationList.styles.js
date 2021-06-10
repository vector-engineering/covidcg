import styled from 'styled-components';

export const MutationListContainer = styled.div`
  display: flex;
  flex-direction: column;
`;

export const MutationListHeader = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
`;

export const OptionSelectContainer = styled.div`
  margin-right: 8px;
  font-weight: normal;
  select {
    margin-left: 0.65em;
    padding: 1px 4px;
    border-radius: 3px;
  }
`;

export const OptionInputContainer = styled.div`
  margin-right: 8px;
  font-weight: normal;
  input {
    max-width: ${({ maxWidth }) => maxWidth};
    margin-left: 0.65em;
    padding: 1px 4px;
    border-radius: 3px;

    border: 1px solid ${({ invalid }) => (invalid ? '#dc3545' : '#aaa')};
    &:focus {
      border: 1px solid ${({ invalid }) => (invalid ? '#dc3545' : '#aaa')};
    }
  }
`;
OptionInputContainer.defaultProps = {
  maxWidth: '4em',
  invalid: false,
};

export const MutationListHeaderTable = styled.table`
  display: grid;
  grid-template-columns: 80px 120px repeat(${({ ncols }) => ncols}, 50px);
  grid-template-rows: 60px auto;
  grid-column-gap: 2px;
  tbody,
  thead,
  tr {
    display: contents;
  }
  min-width: 100%;
  border-bottom: 1px solid #ccc;
`;
MutationListHeaderTable.defaultProps = {
  ncols: 1,
};

export const MutationListHeaderEmpty = styled.th`
  grid-column: span ${({ colSpan }) => colSpan};
`;
MutationListHeaderEmpty.defaultProps = {
  colSpan: 1,
};
export const MutationListHeaderCell = styled.th`
  white-space: nowrap;
  div {
    width: 60px;
    transform: rotate(-45deg) translate(-20px, 20px);
    font-weight: bold;
    text-align: left;
    span {
      padding: 5px 10px;
    }
  }
`;

export const MutationListTable = styled.table`
  display: grid;
  grid-template-columns: 80px 120px repeat(${({ ncols }) => ncols}, 50px);
  grid-column-gap: 2px;
  grid-row-gap: 2px;
  tbody,
  thead,
  tr {
    display: contents;
  }
  min-width: 100%;
  font-weight: normal;
  max-height: 500px;
  overflow-y: scroll;
`;
MutationListTable.defaultProps = {
  ncols: 1,
};

export const MutationRowContainer = styled.tr`
  display: contents;
`;

export const MutationRowBar = styled.td`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-start;

  border-right: 10px solid ${({ barColor }) => barColor};

  grid-row: span ${({ rowSpan }) => rowSpan};
`;
MutationRowBar.defaultProps = {
  rowSpan: 1,
  barColor: '#AAA',
};

export const MutationRowName = styled.td`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-end;
  padding-right: 5px;
  font-family: monospace;
  grid-col: 2;
`;
export const MutationRowHeatmapCellContainer = styled.td`
  padding: 2px 10px;
  text-align: center;
  background-color: ${({ bgColor }) => bgColor};
  font-size: 0.8em;
`;
MutationRowHeatmapCellContainer.defaultProps = {
  bgColor: '#FFF',
};

export const MutationRowHeatmapEmptyCell = styled.td`
  grid-column: span ${({ colSpan }) => colSpan};
`;
MutationRowHeatmapEmptyCell.defaultProps = {
  colSpan: 1,
};
