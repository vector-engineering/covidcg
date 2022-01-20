import styled from 'styled-components';

export const TabIndicatorBG = styled.a`
  display: inline-flex;
  flex-direction: column;
  align-items: stretch;
  width: 130px;
  text-decoration: none;
  color: #444;

  background-color: #eee;
  border: 1px solid #ccc;
  margin-bottom: 0.2em;

  &:hover {
    span {
      background-color: #f8f8f8;
    }
  }
`;

export const TabIndicatorFG = styled.span`
  font-size: 0.85em;
  font-weight: 500;
  text-align: center;
  transition: 0.1s all ease-in-out;
`;
