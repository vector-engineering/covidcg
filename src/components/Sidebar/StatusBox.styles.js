import styled from 'styled-components';

export const Container = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;
  overflow-y: scroll;

  padding: 5px 10px;
  margin-top: 1rem;
  border-top: 1px solid #ccc;
  border-bottom: 1px solid #ccc;
  background-color: #fff;
`;

export const StatusText = styled.div`
  font-size: 0.8rem;
`;

export const Line = styled.p`
  font-size: 1em;
  font-weight: normal;
  margin: 0px;
`;

export const Sequence = styled.span`
  font-family: monospace;
  display: inline;
  margin: 0px;
`;
