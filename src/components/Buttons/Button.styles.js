import styled from 'styled-components';

export const StyledButton = styled.button`
  padding: 5px 10px;
  background-color: ${({ disabled }) => (disabled ? '#ccc' : '#28a745')};
  background-image: ${({ disabled }) =>
    disabled ? 'none' : 'linear-gradient(-180deg, #34d058, #28a745 90%)'};
  color: #fff;
  font-size: 12px;
  border: 1px solid rgba(27, 31, 35, 0.2);
  border-radius: 0.25em;
  outline: none;
  &:active {
    background-color: #279f43;
    background-image: none;
    border-color: rgba(27, 31, 35, 0.5);
  }
  ${({ sticky }) =>
    sticky &&
    `
    &:focus {
      background-color: #279f43;
      background-image: none;
      border-color: rgba(27, 31, 35, 0.5);
    }
  `}
`;
