import styled from 'styled-components';
import { lighten } from 'polished';

import Button from '../Buttons/Button';

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: flex-start;

  width: 600px;
  height: 200px;
`;

export const Content = styled.div`
  width: 100%;
  height: 100%;

  font-size: 1em;
  font-weight: normal;

  display: flex;
  flex-direction: column;
  align-items: stretch;
`;

export const Row = styled.div`
  display: flex;
  flex-direction: row;
  flex-wrap: wrap;

  margin-bottom: 10px;
`;

export const Info = styled.p`
  line-height: normal;
  margin: 5px 0px;
`;

export const RefreshButton = styled(Button)`
  background-image: none;
  font-size: 1rem;
  font-weight: normal;

  background-color: #dc3545;
  color: #fff;

  transition: 0.1s all ease-in-out;

  &:hover {
    background-color: ${() => lighten(0.1, '#dc3545')};
  }
`;
