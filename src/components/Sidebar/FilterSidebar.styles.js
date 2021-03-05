import styled from 'styled-components';

import Button from '../Buttons/Button';

export const Container = styled.div`
  width: 430px;
  display: flex;
  height: 100vh;
  overflow-y: hidden;
`;

export const FilterSidebarContainer = styled.div`
  width: 249px;

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

export const SelectSequencesButton = styled(Button)`
  margin: 5px;
  font-size: 1rem;
`;

export const LegendSidebarContainer = styled.div`
  width: 180px;
  height: 100%;
  border-right: 1px solid #aaa;
  padding-bottom: 15px;
`;
