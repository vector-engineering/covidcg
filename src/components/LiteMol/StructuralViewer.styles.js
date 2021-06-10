import styled from 'styled-components';

export const StructuralViewerContainer = styled.div`
  min-width: 100%;
`;

export const StructuralViewerHeader = styled.div`
  display: flex;
  flex-direction: row;
`;

export const OptionSelectContainer = styled.div`
  margin-right: 8px;
  font-weight: normal;
  select {
    margin-left: 0.65em;
    padding: 1px 4px;
    border-radius: 3px;
  }
`;

export const OptionInputContainer = styled.div`
  margin-right: 8px;
  font-weight: normal;
  input {
    max-width: ${({ maxWidth }) => maxWidth};
    margin-left: 0.65em;
    padding: 1px 4px;
    border-radius: 3px;

    border: 1px solid ${({ invalid }) => (invalid ? '#dc3545' : '#aaa')};
    &:focus {
      border: 1px solid ${({ invalid }) => (invalid ? '#dc3545' : '#aaa')};
    }
  }
`;
OptionInputContainer.defaultProps = {
  maxWidth: '4em',
  invalid: false,
};

export const LiteMolContainer = styled.div`
  min-width: 100%;
  min-height: 600px;
`;
