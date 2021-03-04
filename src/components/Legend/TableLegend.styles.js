import styled from 'styled-components';

export const StyledContainer = styled.div`
  width: 100%;
  height: 100%;
`;

export const Header = styled.div`
  height: 48px;
  border-bottom: 1px solid #ccc;
`;

export const HeaderRow = styled.div`
  width: 100%;
  height: 50%;
  display: flex;
  flex-direction: row;
  align-items: center;
`;

export const StyledColumnHeader = styled.div`
  cursor: pointer;
  width: ${({ width }) => width};
  font-size: 12px;
  padding: 0px 3px;
`;
