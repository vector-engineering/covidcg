import styled from 'styled-components';

const StyledLi = styled.li`
  float: left;
  height: 100%;
  padding: 10px;
`;

export const VOCListContainer = styled.div`
  display: flex;
  flex-direction: column;
  padding: 10px;
`;

export const VOCListHeader = styled.div`
  text-align: center;
`;

export const VOCListTitle = styled.div`
  font-size: 18px;
`;

export const VOCItemTable = styled.table`
  display: table;
`;

export const VOCItemListTitle = styled.span`
  display: inline-block;
  font-size: 16px;
  text-align: center;
  margin-left: auto;
  margin-right: auto;
`;

export const VOCItemList = styled(StyledLi)`
  display: inline;
  float: left;
`;

export const VOCItemContainer = styled.div``;

export const VOCItemName = styled.span``;
export const VOCItemAlias = styled.span``;

export const Td = styled.td`
  vertical-align: top;
`;
