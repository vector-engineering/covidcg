import styled from 'styled-components';

export const OrgLegendContainer = styled.div`
  display: flex;
  flex-direction: column;
  text-align: center;
  margin-left: 10px;
`;

export const OrgLegendTitle = styled.span`
  font-size: 1em;
  font-weight: bold;
`;

export const OrgList = styled.div``;

export const OrgItemContainer = styled.div`
  display: flex;
  flex-direction: column;
`;

export const OrgItem = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  background-color: #dddddd;
  margin: 5px;
  padding: 5px;
`;

export const OrgBadge = styled.div`
  width: 15px;
  height: 15px;
  background-color: ${({ color }) => color};
  border: 2px solid black;
`;

export const OrgName = styled.div`
  margin: auto;
  color: black;
`;
