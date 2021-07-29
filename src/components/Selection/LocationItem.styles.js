import styled from 'styled-components';
import Button from '../Buttons/Button';

export const LocationItemContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  border: 1px solid #ccc;
  border-radius: 3px;
  padding: 1px 5px;
  margin-right: 5px;
`;

export const LocationItemLabel = styled.span`
  margin-right: 3px;
`;

export const DeselectButton = styled(Button)`
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
