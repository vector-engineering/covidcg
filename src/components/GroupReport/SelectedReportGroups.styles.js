import styled from 'styled-components';

export const SelectedReportGroupsContainer = styled.div`
  margin-bottom: 4px;
  font-weight: normal;
`;

export const SelectedReportGroupTitle = styled.div`
  font-size: 1em;
  font-weight: bold;
`;

export const SelectedReportGroupsList = styled.ul`
  margin: 0px;
  padding-left: 5px;
  list-style: none;
`;

export const SelectedReportGroupItemContainer = styled.li`
  display: flex;
  flex-direction: row;
  align-items: center;
  margin-bottom: 2px;
  font-size: 1em;
`;

export const SelectedReportGroupItemTitle = styled.span``;

export const NoReportGroupsSelectedContainer = styled.div``;

export const SelectedReportGroupsButton = styled.button`
  background: none;
  border-style: none;
  color: #aaa;
  cursor: pointer;

  &:hover,
  &:focus {
    color: #ff5555;
  }

  transition: 0.1s all ease-in-out;
`;
