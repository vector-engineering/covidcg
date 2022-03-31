import styled from 'styled-components';

export const WalkthroughListContainer = styled.div`
  padding: 0px 20px;
`;

export const WalkthroughListTitle = styled.div`
  font-size: 1.25rem;
`;

export const WalkthroughItemList = styled.div`
  display: flex;
  flex-direction: row;
`;

export const WalkthroughItemContainer = styled.a`
  min-width: 150px;
  max-width: 300px;
  display: flex;
  flex-direction: column;
  align-items: stretch;
  margin: 0.5rem;
  border: 1px solid #aaa;
  cursor: pointer;
  text-decoration: none;

  transition: 0.1s all ease-in-out;

  &:hover,
  &:active {
    background-color: #f8f8f8;
    text-decoration: underline;
  }
`;

export const WalkthroughItemImage = styled.div`
  height: 2rem;
  background-image: url('${({ image }) => image}');
  background-size: cover;
  background-repeat: no-repeat;
  border-bottom: 1px solid #aaa;
`;

export const WalkthroughItemText = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;
  font-weight: normal;
  padding: 0.25rem;
`;

export const WalkthroughItemTitle = styled.div`
  font-size: 1rem;
  font-weight: 500;
`;

export const WalkthroughItemDescription = styled.div`
  color: #888;
  font-size: 0.8rem;
  text-decoration: none;
`;
