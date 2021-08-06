import styled from 'styled-components';

export const LegendContainer = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;
`;
export const LegendTitle = styled.div``;
export const LegendGradientBar = styled.div`
  width: 100%;
  height: 10px;
  background: ${({ gradient }) => gradient};
`;
export const LegendAxis = styled.div`
  position: relative;
  margin-top: 2px;
  border-top: 1px solid #888;
`;
export const TickContainer = styled.div`
  display: block;
  position: absolute;
  left: ${({ loc }) => loc * 100}%;
`;
export const TickLine = styled.div`
  height: 3px;
  border-left: 1px solid #888;
`;
export const TickLabel = styled.div`
  margin-top: -2px;
  font-size: 0.85em;
  white-space: nowrap;
`;
