import styled from 'styled-components';

const StyledLi = styled.li`
  float: left;
  height: 100%;
  padding: 10px;
`

const VOCGridTitle = styled.span`
  font-size: 16px;
  grid-row: 1;
  justify-self: center;
`

const GridItem = styled(StyledLi)`
  grid-row: auto;
  list-style-type: none;

`

export const VOCGridItem = styled(GridItem)`
  grid-column: 2;
`

export const VOIGridItem = styled(GridItem)`
  grid-column: 5;
`

export const VOCListContainer = styled.div`
  display: flex;
  flex-direction: column;
  padding: 10px;
`

export const VOCListHeader = styled.div`
  text-align: center;
`

export const VOCListTitle = styled.div`
  font-size: 18px;
`

export const VOCItemGrid = styled.div`
  display: grid;
  grid-template-columns: 20px 160px 20px 20px 160px 20px;
  grid-template-rows: auto;
`

export const VOCTitle = styled(VOCGridTitle)`
  grid-column: 2;
`

export const VOITitle = styled(VOCGridTitle)`
  grid-column: 5;
`

export const VOCItemContainer = styled.div`
  &:hover {
    cursor: pointer;
  }
`

export const VOCItemName = styled.span``
export const VOCItemAlias = styled.span``
