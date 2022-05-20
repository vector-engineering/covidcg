import styled from 'styled-components';
import { lighten } from 'polished';

import Button from '../Buttons/Button';

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: flex-start;

  width: 600px;
  height: 400px;
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

export const RadioForm = styled.form`
  display: flex;
  flex-direction: column;
  margin-bottom: 10px;
`;

export const Radio = styled.label`
  display: flex;
  flex-direction: row;
  align-items: center;

  input {
    margin: 0px 8px;
  }
`;

export const CheckboxForm = styled.form`
  display: flex;
  flex-direction: column;
  margin-right: 30px;
`;

export const FormTitle = styled.span`
  font-weight: bold;
  margin-bottom: 10px;
  font-size: 1rem;
`;

export const Checkbox = styled.label`
  display: flex;
  flex-direction: row;
  align-items: center;

  input {
    margin: 0px 8px;
  }

  margin-bottom: 3px;
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

export const SelectInput = styled.label`
  select {
    margin-left: 4px;
  }
`;
