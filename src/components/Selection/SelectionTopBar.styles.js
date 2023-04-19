import styled from 'styled-components';

import Button from '../Buttons/Button';

export const Container = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  padding: 0.25rem 0.5rem;

  grid-column: 2/-1;
  background-color: #f8f8f8;
  border-bottom: 1px solid #ccc;
  font-weight: normal;

  .spacer {
    flex-grow: 1;
  }
`;

export const SelectSequencesButton = styled(Button)`
  font-size: 1rem;
  margin-right: 1rem;
  flex-shrink: 0;
`;

export const StatusBox = styled.div`
  flex-shrink: 1;
`;

export const StatusBlock = styled.div`
  font-size: 0.8rem;
`;
