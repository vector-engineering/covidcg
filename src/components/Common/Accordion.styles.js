import styled from 'styled-components';

export const AccordionTitle = styled.span`
  .question-button {
    font-family: monospace;
    font-size: 1em;
    line-height: normal;

    margin-left: 8px;
    padding: 2px 5px;
    background-color: #fff;
    border: 1px solid #ccc;
    border-radius: 2px;
    &:hover {
      background-color: #f8f8f8;
    }
  }
`;

export const ContentAndButtonWrapper = styled.div`
  display: flex;
  flex-direction: row;
  align-items: flex-start;
  margin-top: 0px;
  margin-bottom: 10px;
  margin-left: 24px;
  margin-right: 24px;
`;

export const ContentWrapper = styled.div`
  width: 100%;
  transition: max-height 0.3s cubic-bezier(0.22, 1, 0.36, 1);
  align-self: center;
  border-radius: 2px;
  max-height: ${({ maxHeight, collapsed }) => {
    if (collapsed) {
      return '1px';
    }
    return maxHeight ? `${maxHeight}` : '100px';
  }};
  overflow-y: scroll;
`;

export const CollapseButton = styled.button`
  border: none;
  cursor: pointer;
  color: #ccc;
  transition: color 0.2s cubic-bezier(0.22, 1, 0.36, 1);
  background-color: rgba(0, 0, 0, 0);
  outline: none;
  width: 32px;
  font-size: 20px;

  &:hover {
    color: black;
  }
`;

export const Title = styled.div`
  padding-left: 12px;
  display: flex;
  align-items: center;
`;

export const HelpButton = styled.button`
  font-size: 0.8em;
  line-height: normal;
  margin-left: 10px;
`;

export const HelpText = styled.div`
  display: ${({ show }) => (show ? 'block' : 'none')};

  margin-bottom: 5px;

  font-weight: normal;
  font-size: 0.9em;
  color: #666;
  line-height: normal;
  p {
    margin-top: 0px;
    margin-bottom: 3px;
  }
`;
HelpText.defaultProps = {
  show: false,
};
