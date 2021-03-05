import styled from 'styled-components';

export const HeaderDiv = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;
  justify-content: flex-start;
  border-bottom: 1px solid #aaa;
  flex-shrink: 0;
`;
export const TitleContainer = styled.div`
  display: flex;
  flex-direction: column;
  align-items: center;
  border-bottom: 1px solid #aaa;
  padding: 5px 0px;

  background-color: #fff;

  h1 {
    font-weight: 700;
    font-size: 1.1em;
    margin: 0px;
    line-height: normal;
  }
`;

export const ImageContainer = styled.div`
  margin-bottom: 2px;
  img {
    width: auto;
    height: 60px;
  }
`;

export const GISAIDContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: center;

  font-size: 0.8rem;
  line-height: normal;

  margin: 5px 12px;

  a {
    display: flex;
    flex-direction: row;
    align-items: center;

    margin-left: 5px;
    img {
      height: 30px;
    }
  }
`;

export const NCBIContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-start;

  padding: 5px 10px;

  span {
    margin-right: 0.5em;
  }

  a {
    display: flex;
    flex-direction: row;
    align-items: center;

    margin-right: 10px;

    img {
      height: 30px;
    }
  }
`;
