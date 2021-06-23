import styled from 'styled-components';

export const StyledButton = styled.button`
  background-color: ${({ disabled }) => (disabled ? '#ccc' : '#2d62fd')};
  background-image: ${({ disabled }) =>
    disabled ? 'none' : 'linear-gradient(-180deg, #22c3ba, #2d62fd 90%)'};
  color: #fff;
  font-size: 12px;
  border: 1px solid rgba(27, 31, 35, 0.2);
  border-radius: 0.25em;
  outline: none;
  margin-left: 10px;

  vertical-align: baseline;

  color: '#ffffff';

  &:active {
    background-color: #2d62fd;
    background-image: none;
    border-color: rgba(27, 31, 35, 0.5);
  }

  ${({ sticky }) =>
    sticky &&
    `
    &:focus {
      background-color: #d5ded6;
      background-image: none;
      border-color: rgba(27, 31, 35, 0.5);
    }
  `}
`;
