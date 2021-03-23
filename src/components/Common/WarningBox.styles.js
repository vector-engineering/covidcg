import styled from 'styled-components';

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
`;
WarningContainer.defaultProps = {
  show: true,
};

export const WarningHeader = styled.div`
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
`;

export const WarningText = styled.p`
  margin: 0px;
  font-weight: normal;
`;
