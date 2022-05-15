import styled from 'styled-components';

import ExternalLink from '../Common/ExternalLink';

export const SelectContainer = styled.div`
  display: flex;
  flex-direction: ${({ direction }) => direction};
  align-items: stretch;
  justify-content: flex-start;

  margin: ${({ direction }) => (direction === 'row' ? '0px' : '5px 0px')};
  padding-right: 10px;
`;

SelectContainer.defaultProps = {
  direction: 'column',
};

export const RadioForm = styled.form`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  margin-right: ${({ direction }) => (direction === 'row' ? '0.75rem' : '0px')};
  margin-bottom: ${({ direction }) => (direction === 'row' ? '0px' : '5px')};

  span.form-title {
    font-size: ${({ direction }) => (direction === 'row' ? 'inherit' : '1rem')};
    font-weight: 500;
    flex-shrink: 0;
  }

  .radio-row {
    display: flex;
    flex-direction: row;
    align-items: center;
    justify-content: flex-start;

    margin-top: ${({ direction }) => (direction === 'row' ? '1px' : '3px')};
    margin-left: ${({ direction }) => (direction === 'row' ? '0px' : '12px')};
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

RadioForm.defaultProps = {
  direction: 'column',
};

export const Link = styled(ExternalLink)`
  font-size: 0.9em;
  margin-left: 10px;
`;

export const HintText = styled.div`
  font-size: 1em;
  color: #888;
`;

export const ReferenceSelectRow = styled.div`
  display: flex;
  flex-direction: column;
  align-items: flex-start;

  select {
    margin-top: 2px;
    width: 100%;
  }

  .reference-description {
    font-size: 0.8em;
    font-weight: normal;
  }
`;
