import styled from 'styled-components';

export const PlotHeader = styled.div`
  display: grid;
  grid-template-columns: minmax(150px, 1fr) 3fr;
  grid-template-rows: repeat(3, auto);
  grid-row-gap: 3px;
`;

export const PlotOptionsRow = styled.div`
  grid-column: 2;
  display: flex;
  flex-direction: row;
  align-items: center;
`;
