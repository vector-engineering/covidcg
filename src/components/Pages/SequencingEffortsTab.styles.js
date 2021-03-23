import styled from 'styled-components';

export const Container = styled.div`
  padding-top: 20px;
  overflow-x: hidden;
`;

export const Header = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  width: 100%;

  padding: 0px 20px;

  p {
    font-weight: normal;
    line-height: normal;
    max-width: 800px;
    margin: 5px 0px;
  }
`;
export const Title = styled.h2`
  font-size: 1.5em;
  font-weight: bold;

  margin: 0px;
  margin-bottom: 10px;
`;
