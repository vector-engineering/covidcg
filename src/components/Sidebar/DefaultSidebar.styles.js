import styled from 'styled-components';

export const SidebarContainer = styled.div`
  position: fixed;
  top: 0;
  width: 249px;

  background-color: #f8f8f8;
  //padding-right: 10px;
  padding-bottom: 15px;
  border-right: 1px solid #aaa;
  display: flex;
  flex-direction: column;
  height: 100vh;
  overflow-y: hidden;

  .filter-sidebar-tooltip {
    background-color: #fff;
    font-weight: normal;
    p {
      margin-top: 2px;
      margin-bottom: 2px;
    }
  }
`;

export const SidebarChunk = styled.div`
  margin: 5px 12px;
  padding: 5px 0px;
  font-weight: normal;
`;
