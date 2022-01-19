import styled from 'styled-components';

export const ContainerDiv = styled.div`
  width: 100%;
  height: 100%;

  margin-top: 2px;
  padding-top: 8px;

  display: flex;
  flex-direction: column;
  // overflow-y: hidden;
  overflow-y: scroll;
`;

export const DropdownHeader = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;

  margin-left: 15px;
`;

export const Title = styled.span`
  font-size: 1rem;
  font-weight: 500;
`;

export const UnselectButton = styled.button`
  display: ${({ show }) => (show ? 'block' : 'none')};
  margin-left: 10px;
`;
UnselectButton.defaultProps = {
  show: true,
};

export const SelectedLocationsContainer = styled.div`
  display: flex;
  flex-direction: row;
  flex-wrap: wrap;
  padding: 5px 0px;
`;

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
