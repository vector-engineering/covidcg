import styled from 'styled-components';

export const HomePageDiv = styled.div`
  display: grid;
  grid-template-columns: [col1] 250px [col2] 180px [col3] auto [col4];
  grid-template-rows: [row1] 50px [row2] auto [row3];
  width: 100vw;
  height: 100vh;
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

export const LegendContainer = styled.div`
  grid-row: 2/3;
  width: 180px;
  height: 100%;
  border-right: 1px solid #aaa;
  padding-bottom: 15px;
`;

export const PlotContainer = styled.div`
  grid-column: ${({ showDefaultSidebar }) =>
    showDefaultSidebar ? '2/-1' : '3/-1'};
  grid-row: ${({ showDefaultSidebar }) =>
    showDefaultSidebar ? '1/-1' : '2/-1'};
  display: flex;
  flex-direction: column;
  width: 100%;
  height: 100vh;
  box-sizing: border-box;
  position: relative;
  overflow-y: scroll;
  border-left: 1px #eaeaea solid;
`;
