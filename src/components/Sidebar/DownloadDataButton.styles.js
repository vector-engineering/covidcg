import styled from 'styled-components';

import Button from '../Buttons/Button';

export const ButtonContainer = styled.div`
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  margin-bottom: 5px;
`;

export const DownloadButton = styled(Button)`
  background-color: ${({ disabled }) => (disabled ? '#ccc' : '#eee')};
  background-image: none;
  color: #000;
  border-color: #666;

  .caret:after {
    border-top-color: #666;
  }

  &:hover,
  &:focus {
    background-color: ${({ disabled }) => (disabled ? '#ccc' : '#ddd')};
  }
`;
