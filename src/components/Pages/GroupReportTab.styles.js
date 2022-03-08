import styled from 'styled-components';

const MEDIAWIDTH = '1475px';

export const GroupReportTabContainer = styled.div`
  display: grid;
  grid-template-rows: auto;
  grid-template-columns: auto
  max-height: 100vh;

  @media (max-width: ${MEDIAWIDTH}) {
    grid-template-rows: auto;
    grid-template-columns: 250px minmax(100px, 1fr);
    overflow-y: auto;
  }
`;

export const GroupTreePlotContainer = styled.div`
  grid-row: 1 / 3;
  grid-column: 1;

  @media (max-width: ${MEDIAWIDTH}) {
    grid-row: 1 / 4;
  }
`;

export const MutationsContainer = styled.div`
  grid-column: 2;
  grid-row: 1 / 3;

  @media (max-width: ${MEDIAWIDTH}) {
    grid-row: 2;
  }
`;

export const VOCContainer = styled.div`
  grid-column: 3;
  grid-row: 1;

  @media (max-width: ${MEDIAWIDTH}) {
    grid-column: 2;
    grid-row: 1;
  }
`;

export const StructuralViewerContainer = styled.div`
  grid-column: 3;
  grid-row: 2;

  @media (max-width: ${MEDIAWIDTH}) {
    grid-column: 2;
    grid-row: 3;
  }
`;
