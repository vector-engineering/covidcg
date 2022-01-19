import styled from 'styled-components';

export const GroupSelectContainer = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  padding-left: 15px;
  margin-bottom: 10px;

  span.title {
    font-weight: 500;
    font-size: 1rem;
    margin-bottom: 5px;
  }
`;

export const GroupSelectForm = styled.div``;

export const GroupSelectedItemsContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: flex-start;
`;

export const GroupTitle = styled.span`
  font-weight: 500;
  margin-right: 0.5rem;
`;

export const GroupSelectedItemList = styled.div`
  display: flex;
  flex-direction: row;
  flex-wrap: wrap;
  align-items: flex-start;
`;

export const GroupSelectedItemContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  border: 1px solid #ccc;
  border-radius: 3px;
  padding: 1px 5px;
  margin-right: 5px;
`;

export const GroupSelectedItemLabel = styled.span`
  margin-right: 3px;
`;
