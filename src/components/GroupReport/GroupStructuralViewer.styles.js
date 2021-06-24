import styled from 'styled-components';
import Button from '../Buttons/Button';

export const StructuralViewerContainer = styled.div`
  min-width: 100%;
`;

export const StructuralViewerHeader = styled.div`
  display: flex;
  flex-direction: row;

  .spacer {
    flex-grow: 1;
  }
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

export const InvalidText = styled.span`
  margin: 0px 5px;
  font-size: 0.9em;
  font-weight: normal;
  line-height: normal;
  color: #dc3545;
`;
export const ConfirmButton = styled(Button)`
  margin-left: 5px;
`;

export const LiteMolContainer = styled.div`
  min-width: 100%;
  min-height: 600px;
  margin-top: 5px;
`;
