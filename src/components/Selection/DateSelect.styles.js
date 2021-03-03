import styled from 'styled-components';

export const Container = styled.div`
  display: flex;
  flex-direction: column;

  padding-left: 15px;
  margin-top: 10px;
  margin-bottom: 10px;

  span.title {
    font-weight: 500;
    font-size: 1rem;
    margin-bottom: 5px;
  }
`;

export const PresetForm = styled.div`
  margin: 5px 0px;

  label {
    margin-right: 5px;
  }

  select {
  }
`;

export const DateForm = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
`;

export const FormColumn = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  margin-right: 10px;

  input {
    font-family: inherit;
  }
`;
