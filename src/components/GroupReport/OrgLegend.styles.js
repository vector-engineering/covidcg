import styled from 'styled-components';

export const OrgLegendContainer = styled.div`
  display: flex;
  flex-direction: row;
  text-align: center;
  margin-left: 10px;
  margin-bottom: 5px;
  justify-content: center;
  align-items: center;
`;

export const OrgLegendTitle = styled.span`
  font-size: 1em;
  font-weight: bold;
  margin-right: 5px;
`;

export const OrgList = styled.div``;

export const OrgItemContainer = styled.div`
  display: flex;
  flex-direction: row;
`;

export const OrgItem = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  margin: auto;
  padding: 5px;
`;

export const OrgBadge = styled.div`
  width: 15px;
  height: 15px;
  background-color: ${({ color }) => color};
`;

export const OrgName = styled.div`
  margin-left: 5px;
  color: black;
`;
