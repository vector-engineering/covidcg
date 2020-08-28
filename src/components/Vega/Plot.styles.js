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

export const OptionInputContainer = styled.div`
  margin-right: 8px;
  font-weight: normal;
  input {
    max-width: ${({ maxWidth }) => maxWidth};
    margin-left: 0.65em;
    padding: 1px 4px;
    border-radius: 3px;
  }
`;
OptionInputContainer.defaultProps = {
  maxWidth: '4em',
};

export const OptionCheckboxContainer = styled.div`
  label {
    display: flex;
    flex-direction: row;
    align-items: center;

    font-weight: normal;
    input {
      margin: 0px 6px 0px 0px;
      padding: 0px;
      border-radius: 3px;
    }
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

export const WarningContainer = styled.div`
  display: ${({ show }) => (show ? 'flex' : 'none')};
  flex-direction: column;
  align-items: stretch;
  justify-content: flex-start;

  // colors from Bootstrap
  background-color: #fff3cd;
  border: 1px solid #aaa;
  border-radius: 5px;

  margin: 0px 10px;
  margin-bottom: 15px;
  padding: 10px 20px;

  .warning-header {
    display: flex;
    flex-direction: row;
    align-items: center;
    margin-bottom: 5px;
    .warning-title {
      font-size: 1.25em;
    }
    .spacer {
      flex-grow: 1;
    }
  }

  .warning-text {
    margin: 0px;
    font-weight: normal;
  }
`;
WarningContainer.defaultProps = {
  show: true,
};
