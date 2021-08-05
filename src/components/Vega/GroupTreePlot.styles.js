import styled from 'styled-components';

export const TreePlotContainer = styled.div`
  width: ${({ width }) => width}px;
  height: 100vh;

  display: flex;
  flex-direction: column;
  align-items: stretch;

  position: sticky;
  left: 0px;
  top: 0px;

  font-weight: normal;
`;
TreePlotContainer.defaultProps = {
  width: 300,
};

export const Header = styled.div`
  padding: 5px;
  height: ${({ headerHeight }) => headerHeight}px;
`;

export const HeaderRow = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
`;

export const Title = styled.span`
  font-size: 1.2em;
  font-weight: bold;
`;
export const SubTitle = styled.div`
  font-size: 0.8em;
  line-height: normal;
  margin-bottom: 6px;
`;

export const SelectContainer = styled.div`
  select {
    margin-left: 0.5em;
  }
`;

export const TreeScrollContainer = styled.div`
  overflow-y: scroll;
  height: calc(100vh - ${({ headerHeight }) => headerHeight}px);
`;
