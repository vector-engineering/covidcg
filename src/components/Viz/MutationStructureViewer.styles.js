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
