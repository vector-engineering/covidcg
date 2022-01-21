import styled from 'styled-components';

export const Name = styled.span`
  display: inline-flex;
  flex-wrap: wrap;
  align-items: center;
  font-size: 1em;
  height: 100%;
  color: ${({ selected }) => (selected ? 'white' : 'black')};
  pointer-events: none;

  &:after {
    margin-left: 5px;
    font-family: monospace;
    content: ${({ whoLabel }) => (whoLabel ? JSON.stringify(whoLabel) : '')};
  }
`;
