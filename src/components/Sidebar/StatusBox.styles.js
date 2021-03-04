import styled from 'styled-components';

import Button from '../Buttons/Button';

export const Container = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;
  overflow-y: scroll;

  padding: 5px 10px;
  border-top: 1px solid #ccc;
  border-bottom: 1px solid #ccc;
  background-color: #fff;
`;

export const StatusText = styled.div`
  font-size: 0.8rem;
`;

export const ButtonContainer = styled.div`
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  margin-top: 5px;
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

export const DownloadButton = styled(Button)`
  background-color: #eee;
  background-image: none;
  color: #000;
  border-color: #666;

  .caret:after {
    border-top-color: #666;
  }

  &:hover,
  &:focus {
    background-color: #ddd;
  }
`;
