import styled from 'styled-components';
import Button from './Button';

const DeselectButton = styled(Button)`
  background-color: #fff;
  background-image: none;
  color: #888;
  padding: 0px 5px;
  border: none;
  cursor: pointer;
  transition: 0.1s all ease-in-out;

  &:hover,
  &:active {
    background-color: #f8f8f8;
    color: #ff5555;
  }
`;

export default DeselectButton;
