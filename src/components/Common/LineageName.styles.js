import styled from 'styled-components';

export const Name = styled.span`
  margin-left: 5px;
  color: ${({ selected }) => (selected ? 'white' : 'black')};

  &:after {
    margin-left: 5px;
    font-family: monospace;
    content: ${({ whoLabel }) => (whoLabel ? JSON.stringify(whoLabel) : '')};
  }
`;
