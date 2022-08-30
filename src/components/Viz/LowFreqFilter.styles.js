import styled from 'styled-components';
import { lighten } from 'polished';
import Button from '../Buttons/Button';

export const LowFreqFilterContainer = styled.div`
  display: grid;
  grid-gap: 5px;
  grid-auto-flow: column;
  align-items: center;
`;

export const LowFreqFilterSelectContainer = styled.label`
  margin-right: 8px;
  font-weight: normal;
  select {
    margin-left: 0.65em;
    padding: 1px 4px;
    border-radius: 3px;
  }
`;

export const LowFreqFilterInputContainer = styled.label`
  margin-right: 8px;
  font-weight: normal;
  input {
    max-width: ${({ maxWidth }) => maxWidth};
    padding: 1px 4px;
    border-radius: 3px;

    border: 1px solid ${({ invalid }) => (invalid ? '#dc3545' : '#aaa')};
    &:focus {
      border: 1px solid ${({ invalid }) => (invalid ? '#dc3545' : '#aaa')};
    }
  }
`;
LowFreqFilterInputContainer.defaultProps = {
  maxWidth: '4em',
  invalid: false,
};

export const ApplyButton = styled(Button)`
  background-image: none;
  font-size: 0.85em;
  font-weight: normal;

  background-color: ${({ invalid }) => (invalid ? '#DDD' : '#28a745')};
  color: ${({ invalid }) => (invalid ? '#888' : '#FFF')};

  transition: 0.1s all ease-in-out;

  &:hover {
    background-color: ${({ invalid }) =>
      invalid ? '#DDD' : lighten(0.1, '#28a745')};
  }
`;
ApplyButton.defaultProps = {
  invalid: false,
};
