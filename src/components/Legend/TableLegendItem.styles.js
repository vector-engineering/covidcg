import styled from 'styled-components';

export const Container = styled.div`
  border-bottom: 1px solid #eee;
  display: flex;
  align-items: center;
  ${({ hovered, selected }) => {
    if (selected) return 'background-color: rgba(0,0,0,0.3);';
    if (hovered) return 'background-color: rgba(0,0,0,0.1);';
  }}
  font-size: 12px;
`;

export const ColorBar = styled.div`
  height: 100%;
  border-right: 5px solid ${({ color }) => color};
  margin-right: 6px;
`;

export const GroupNameContainer = styled.div`
  border-right: 1px #eee solid;
  width: 40%;
  height: 100%;
  display: flex;
  align-items: center;

  font-size: 0.85em;
  line-height: normal;
  word-break: break-all;
`;

export const PercentageContainer = styled.div`
  width: 30%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: flex-start;

  padding: 0px 5px;
  font-size: 0.7rem;
`;
