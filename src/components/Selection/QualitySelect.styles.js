import styled from 'styled-components';

export const QualitySelectContainer = styled.div`
  display: flex;
  flex-direction: column;

  padding-left: 15px;
  margin-bottom: 10px;

  span.title {
    font-weight: 500;
    font-size: 1rem;
    margin-bottom: 5px;
  }
`;

export const FormRow = styled.div`
  display: flex;
  flex-direction: row;
  align-items: flex-end;
`;

export const TitleColumn = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  margin-right: 10px;
  font-weight: 500;
`;

export const FormColumn = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;
  max-width: 6em;

  margin-right: 10px;

  input {
    font-family: inherit;
  }
`;
