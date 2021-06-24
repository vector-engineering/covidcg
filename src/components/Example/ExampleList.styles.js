import styled from 'styled-components';

export const ExampleListContainer = styled.div``;

export const ExampleHeader = styled.div`
  padding-left: 10px;
  max-width: 800px;

  p {
    font-weight: normal;
    line-height: normal;
  }
`;

export const ExampleTitle = styled.h2``;

export const ExampleItemList = styled.div`
  display: flex;
  flex-direction: row;
  flex-wrap: wrap;
`;
export const ExampleItem = styled.a`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  width: 300px;
  height: 250px;
  margin: 10px;
  border: 1px solid #ccc;
  box-shadow: 0px 3px 4px #aaa;

  text-decoration: none;

  transition: 0.1s all ease-in-out;
  &:hover {
    border: 1px solid #00e;
  }
`;
export const ExampleItemImage = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: center;
  overflow: hidden;

  img {
    width: auto;
    height: 150px;
  }
`;

export const ExampleItemFooter = styled.div`
  height: 100px;
  padding: 10px;

  border-top: 1px solid #aaa;

  .example-item-title {
    display: block;
    font-size: 1.1em;
    font-weight: 700;
    color: #444;
  }

  .example-item-description {
    display: block;
    font-size: 0.85em;
    font-weight: normal;
    color: #888;
    margin: 3px 0px;
    line-height: normal;
  }
  //
`;
