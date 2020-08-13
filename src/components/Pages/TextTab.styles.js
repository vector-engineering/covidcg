import styled from 'styled-components';

export const TabContainer = styled.div`
  background-color: #f8f8f8;
`;

export const Content = styled.div`
  flex-grow: 1;
  max-width: 1000px;
  padding: 20px;
  margin: 0 auto;

  background-color: #fff;
`;

export const ContentSection = styled.div`
  border-top: 1px solid #ccc;
  font-weight: normal;
  font-size: 1em;
  padding-top: 10px;
  padding-bottom: 10px;

  .section-title {
    font-weight: bold;
    font-size: 1.5em;
  }

  .content-block {
    display: flex;
    flex-direction: row;
    align-items: flex-start;

    margin-top: 10px;

    .content-text {
      width: 60%;
      padding-right: 10px;

      .content-subtitle {
        display: block;
        font-size: 1.25em;
        font-weight: bold;
        margin-bottom: 10px;
      }

      p {
        margin-top: 5px;
        margin-bottom: 5px;
      }
    }
    .content-images {
      width: 40%;
      display: flex;
      flex-direction: column;
      align-items: flex-start;
      justify-content: flex-start;

      padding-left: 10px;

      div {
        margin-bottom: 5px;
      }
    }
  }
`;

export const ImageRow = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-start;
  a {
    margin-right: 10px;
  }
`;
