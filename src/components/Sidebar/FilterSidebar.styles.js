import styled from 'styled-components';

export const FilterSidebarContainer = styled.div`
  grid-row: 1/-1;
  width: 249px;
  max-height: 100vh;

  background-color: #f8f8f8;
  //padding-right: 10px;
  padding-bottom: 15px;
  border-right: 1px solid #aaa;
  display: flex;
  flex-direction: column;
  overflow-y: hidden;
  height: 100%;

  .filter-sidebar-tooltip {
    background-color: #fff;
    font-weight: normal;
    p {
      margin-top: 2px;
      margin-bottom: 2px;
    }
  }
`;
