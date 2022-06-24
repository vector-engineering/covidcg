import styled from 'styled-components';
import { lighten } from 'polished';

import Button from '../Buttons/Button';

export const Overlay = styled.div`
  visibility: ${({ visible }) => (visible ? 'visible' : 'hidden')};
  position: absolute;
  top: 0;
  left: 0;
  width: 100%;
  height: 100%;
  z-index: 2;
  background-color: rgba(230, 230, 230, 0.8);

  display: flex;
  flex-direction: column;
  justify-content: center;
`;

export const ProgressContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: center;
`;

export const ProgressText = styled.span`
  margin-left: 20px;
  font: 14px 'Helvetica Neue', Helvetica, Arial, sans-serif;
  line-height: 1.4em;
  color: #4d4d4d;

  font-size: 2rem;
  font-weight: 500;
`;

export const HeaderContainer = styled.div`
  position: fixed;
  top: 0;
  left: 0;
  width: 100%;
  flex-shrink: 0;
  background-color: #fff;

  height: 40px;

  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-start;

  border-bottom: 1px solid #ccc;
  box-shadow: 0px 0px 5px 1px #ddd;
`;

export const HeaderRow = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  height: 100%;
  width: 100%;
  padding: 3px 5px;
`;

export const HeaderButtons = styled.div`
  display: flex;
  flex-direction: row;
`;

export const TitleContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: flex-start;
  margin-left: 20px;

  .title {
    h2 {
      margin-bottom: 0px;
      margin-top: 0px;
    }
  }
  .spacer {
    flex-grow: 1;
  }
  .close-button {
  }
`;

export const CancelButton = styled(Button)`
  background-color: #ddd;
  color: #000;
  background-image: none;
  font-size: 1rem;
  font-weight: normal;
  margin-right: 10px;

  &:hover,
  &:active {
    background-color: #eee;
  }
`;

export const InvalidText = styled.span`
  display: flex;
  flex-direction: row;
  align-items: center;
  margin-right: 20px;

  font-size: 1rem;
  font-weight: normal;
  line-height: normal;
  color: #dc3545;
}
`;

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: flex-start;

  width: ${({ width }) => width}px;
  height: ${({ height }) => height}px;
`;
Wrapper.defaultProps = {
  width: 500,
  height: 200,
};

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

export const FormTitle = styled.span`
  font-weight: bold;
  margin-bottom: 10px;
  font-size: 1rem;
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

export const SelectInput = styled.label`
  select {
    margin-left: 4px;
  }
`;

export const TextInput = styled.label`
  input {
    margin-left: 4px;
  }
`;

export const CheckboxForm = styled.form`
  display: flex;
  flex-direction: column;
  margin-right: 30px;
`;

export const CheckboxInput = styled.label`
  input {
    margin-right: 4px;
  }
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
