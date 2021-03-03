import styled from 'styled-components';

export const SelectContainer = styled.div`
  display: flex;
  flex-direction: column;

  padding-left: 15px;

  span.title {
    font-weight: 500;
    font-size: 1rem;
    margin-bottom: 5px;
  }
`;

export const SelectItem = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;

  padding-left: 10px;
  margin-bottom: 3px;
`;

export const Radio = styled.label`
  display: flex;
  flex-direction: row;
  align-items: center;

  flex-shrink: 0;
  margin-right: 5px;
`;

export const RadioLabel = styled.span`
  color: ${({ itemSelected }) => (itemSelected ? 'inherit' : '#AAA')};
  margin-left: 3px;
`;
RadioLabel.defaultProps = {
  itemSelected: false,
};

export const NumberInput = styled.label`
  flex-shrink: 0;
  width: 50px;
  margin-right: 5px;

  input {
    width: 100%;
  }
`;
