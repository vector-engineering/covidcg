import styled from 'styled-components';

export const ListContainer = styled.div`
  margin-bottom: 5px;

  table {
    width: 100%;
    display: grid;
    border: 1px solid #aaa;

    grid-template-columns: 2fr 1fr 6fr 2fr 2fr;
    background-color: #ccc;
    grid-column-gap: 1px;
    grid-row-gap: 1px;
  }

  tbody {
    display: contents;
  }
  tr {
    display: contents;
  }

  td,
  th {
    background-color: #fff;
    text-align: left;
    padding: 2px 5px;
  }
  td {
    font-weight: normal;
  }
`;

export const ItemContainer = styled.tr``;
