import styled from 'styled-components';

export const HomeTabContainer = styled.div`
  background-color: #f8f8f8;
`;

export const HomeTabContent = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  flex-grow: 1;
  padding: 20px;
  margin: 0 auto;

  background-color: #fff;
`;

export const PubBanner = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  background-color: #fff3cd;
  padding: 5px;
  border-bottom: 1px solid #aaa;

  div {
    margin: 0px;
    margin-left: auto;
  }
`;

export const CloseButton = styled.button`
  margin-left: auto;
`;
