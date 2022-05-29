import styled from 'styled-components';

export const StyledContainer = styled.div`
  width: 100%;
  height: 100%;
`;

export const Header = styled.div`
  border-bottom: 1px solid #ccc;
`;

export const HeaderRow = styled.div`
  width: 100%;
  height: 24px;
  display: flex;
  flex-direction: row;
  align-items: center;

  font-size: 12px;
`;

export const StyledColumnHeader = styled.div`
  cursor: pointer;
  width: ${({ width }) => width};
  padding: 0px 3px;
`;
