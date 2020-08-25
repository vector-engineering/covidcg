import styled from 'styled-components';

export const PlotOptions = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-start;
  margin-bottom: 10px;
  padding-right: 10px;

  .spacer {
    flex-grow: 1;
  }
`;

export const OptionSelectContainer = styled.div`
  margin-right: 8px;
  font-weight: normal;
  select {
    margin-left: 0.65em;
    padding: 1px 4px;
    border-radius: 3px;
  }
`;

export const PlotTitle = styled.div`
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  justify-content: center;

  border-right: 1px solid #ccc;
  margin-right: 10px;
  padding-right: 10px;
  padding-left: 18px;

  line-height: normal;

  .title {
    font-size: 1.25em;
  }
  .subtitle {
    font-size: 0.9em;
    font-weight: normal;
  }
`;
