import styled from 'styled-components';

export const HomePageDiv = styled.div`
  display: grid;
  grid-template-columns: [col1] 250px [col2] 180px [col3] auto [col4];
  grid-template-rows: [row1] auto [row2];
  width: 100vw;
  position: relative;
  overflow-y: hidden;

  .main-tooltip {
    max-width: 300px;

    background-color: #fff;
    font-weight: normal;
    p {
      margin-top: 2px;
      margin-bottom: 2px;
    }
  }
`;

export const PlotContainer = styled.div`
  grid-column: ${({ showDefaultSidebar }) =>
    showDefaultSidebar ? 'col2 / col4' : 'col3 / col4'};
  grid-row: row1 / row2;
  display: flex;
  flex-direction: column;
  width: 100%;
  height: 100vh;
  box-sizing: border-box;
  position: relative;
  overflow-y: scroll;
  border-left: 1px #eaeaea solid;
`;
