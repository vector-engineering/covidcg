import styled from 'styled-components';

import ExternalLink from '../Common/ExternalLink';

export const SelectContainer = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;
  justify-content: flex-start;

  margin: 5px 0px;
  padding-right: 10px;
`;

export const RadioForm = styled.form`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  margin-bottom: 5px;

  span.form-title {
    font-size: 1rem;
    font-weight: 500;
    flex-shrink: 0;
  }

  .radio-row {
    display: flex;
    flex-direction: row;
    align-items: center;
    justify-content: flex-start;

    margin-top: 3px;
    margin-left: 12px;
  }

  .radio-item {
    margin-left: 0.65em;
    &:first-child {
      margin-left: 0px;
    }

    display: flex;
    flex-direction: row;
    align-items: center;

    input {
      margin: 0px;
      margin-right: 0.5em;
    }

    span.disabled-text {
      font-weight: normal;
      font-size: 0.8em;
      color: #888;
      flex-shrink: 1;
      white-space: normal;
      line-height: normal;
      margin-left: 5px;
    }
  }
`;

export const Link = styled(ExternalLink)`
  font-size: 0.9em;
  margin-left: 10px;
`;

export const HintText = styled.div`
  font-size: 1em;
  color: #888;
`;
