import styled from 'styled-components';
import Button from '../Buttons/Button';

export const MutationStructureViewerContainer = styled.div``;

export const LiteMolContainer = styled.div`
  min-width: 100%;
  min-height: 500px;
  // margin-top: 5px;
`;

export const InvalidText = styled.span`
  margin: 0px 5px;
  font-size: 0.9em;
  font-weight: normal;
  line-height: normal;
  color: #dc3545;
`;
export const ConfirmButton = styled(Button)`
  margin-left: 2px;
  margin-right: 10px;
`;

export const HighlightedMutations = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;

  ul {
    list-style-type: none;
    display: flex;
    flex-direction: row;
    align-items: center;
    padding-left: 0.5rem;
    font-weight: normal;
    font-family: monospace;
    margin: 0.25rem 0rem;

    li {
      margin-left: 0.25rem;
      margin-right: 0.25rem;
    }
  }
`;
