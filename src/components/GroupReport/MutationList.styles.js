import styled from 'styled-components';
import { transparentize } from 'polished';

export const MutationListContainer = styled.div`
  display: flex;
  flex-direction: column;
  position: sticky;
  left: 0px;
  top: 0px;
  align-items: stretch;

  height: 99vh;

  padding-top: 4px;
`;

export const MutationListHeader = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;

  .spacer {
    flex-grow: 1;
  }
`;

export const MutationListTitle = styled.div`
  font-weight: 500;
`;

export const MutationContentContainer = styled.div``;

export const MutationInnerContainer = styled.div`
  max-width: 100%;
  overflow-x: auto;
`;

export const OptionSelectContainer = styled.div`
  margin-right: 12px;
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

export const OptionCheckboxContainer = styled(OptionInputContainer)`
  input {
    margin-right: 0.4em;
  }
`;

export const MutationListHeaderTable = styled.table`
  position: sticky;
  top: 0px;
  z-index: 1;
  background-color: white;
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

export const DeleteButton = styled.td`
  background: none;
  border-style: none;
  color: #aaa;
  cursor: pointer;
  text-align: center;

  &:hover,
  &:focus {
    color: #ff5555;
  }

  transition: 0.1s all ease-in-out;
`;

export const MutationListTable = styled.table`
  position: relative;
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

  cursor: pointer;
  transition: 0.1s all ease-in-out;

  &:hover,
  &:active {
    border-right-color: ${({ barColor }) => transparentize(0.4, barColor)};
  }
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
  grid-column: 2 / span ${({ colSpan }) => colSpan};
`;
MutationRowName.defaultProps = {
  colSpan: 1,
};

export const MutationRowHeatmapCellContainer = styled.td`
  padding: 2px 10px;
  text-align: center;
  background-color: ${({ bgColor }) => bgColor};
  font-size: 0.8em;
`;
MutationRowHeatmapCellContainer.defaultProps = {
  bgColor: '#FFF',
};
