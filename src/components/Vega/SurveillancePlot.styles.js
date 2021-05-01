import styled from 'styled-components';
import PropTypes from 'prop-types';

export const PlotContainer = styled.div`
  overflow-x: hidden;
`;

export const PlotTitle = styled.div`
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  justify-content: center;

  margin-right: 10px;
  padding-right: 10px;
  padding-left: 18px;

  line-height: normal;

  .title {
    font-size: 1.25rem;
  }
  .subtitle {
    font-size: 0.9em;
    font-weight: normal;
  }
`;

export const HelpText = styled.p`
  margin: 0px 20px;
  margin-bottom: 5px;

  font-weight: normal;
  font-size: 14px;
  line-height: normal;
`;

export const OptionsColumn = styled.div`
  display: flex;
  flex-direction: column;
`;

export const OptionsRow = styled.div`
  display: flex;
  flex-direction: row;
  margin-top: 5px;
`;

export const PlotAndLegend = styled.div`
  display: flex;
  flex-direction: row;
`;

export const Legend = styled.div`
  min-width: 150px;
  display: flex;
  flex-direction: column;
  padding-right: 5px;
`;

export const LegendTitle = styled.span`
  margin-bottom: 3px;
`;

export const LegendList = styled.div`
  display: flex;
  flex-direction: column;
  max-height: 400px;
  overflow-y: scroll;
`;

export const LegendItemContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  font-weight: normal;
  font-size: 0.8rem;
  padding: 2px 5px;
  border-top: 1px solid #eee;
  border-bottom: 1px solid #eee;
  border-right: 1px solid #eee;
  border-left: 6px solid ${({ barColor }) => barColor};
  cursor: default;
  word-break: break-all;
  color: ${({ hovered }) => (hovered ? '#000' : '#666')};
  background-color: ${({ hovered }) => (hovered ? '#F8F8F8' : '#FFF')};
`;
LegendItemContainer.propTypes = {
  hovered: PropTypes.bool,
  barColor: PropTypes.string,
};

export const LegendItemName = styled.div`
  width: 100px;
`;
export const LegendItemCounts = styled.div`
  display: flex;
  flex-direction: row;
  align-items: flex-start;

  width: 30px;
  padding-left: 3px;
  border-left: 1px solid #eee;
`;
