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

export const ApplyButton = styled(Button)`
  background-image: none;
  font-size: 1rem;
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
