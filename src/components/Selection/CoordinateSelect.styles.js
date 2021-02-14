import styled from 'styled-components';

import Button from '../Buttons/Button';
import DropdownTreeSelect from 'react-dropdown-tree-select';

export const SelectContainer = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;
  justify-content: flex-start;

  margin: 5px 0px;
  padding-right: 10px;

  span.title {
    font-weight: 500;
    font-size: 1rem;
    margin-bottom: 5px;
  }
`;

export const ModeSelectForm = styled.div``;

export const ModeRadioHorizontal = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;

  margin-left: 10px;
  margin-bottom: 5px;
`;

export const ModeRadioVertical = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;
  justify-content: flex-start;

  margin-left: 10px;
  margin-bottom: 8px;
`;

export const ModeLabel = styled.label`
  display: flex;
  flex-direction: row;
  align-items: center;

  padding-right: 10px;

  input.radio-input {
    margin: 0px 8px 0px 0px;
  }

  span.select-text {
    flex-shrink: 0;
  }

  span.hint-text {
    font-size: 0.85em;
    color: #888;
    margin-left: 8px;
    font-weight: normal;
    line-height: normal;
  }
`;

export const SelectForm = styled.form`
  display: flex;
  flex-direction: row;
  align-items: center;
  flex-grow: 1;

  margin-left: 8px;

  label {
    display: flex;
    flex-direction: row;
    align-items: center;
    justify-content: flex-start;
  }

  select {
    background-color: white;
    flex-grow: 1;
    padding: 1px 5px;
    width: 100%;
    border-radius: 3px;
  }
`;

export const DomainSelectForm = styled.form`
  display: flex;
  flex-direction: row;
  align-items: center;
  flex-grow: 1;

  margin-top: 3px;
  margin-left: 20px;

  span {
    font-weight: normal;
    margin-right: 8px;
  }

  label {
    display: flex;
    flex-direction: row;
    align-items: center;
    justify-content: flex-start;
  }

  select {
    background-color: white;
    flex-grow: 1;
    padding: 1px 5px;
    width: 100%;
    border-radius: 3px;
  }
`;

export const PrimerSelectContainer = styled.div`
  span.placeholder {
    &:after {
      content: '${(props) => props.placeholderText}';
    }
  }
`;

PrimerSelectContainer.defaultProps = {
  placeholderText: '',
};

export const PrimerSelect = styled(DropdownTreeSelect)`
  a.dropdown-trigger {
    &:focus {
      outline: none;
    }

    ul.tag-list {
      margin-top: 6px;
      margin-bottom: 0px;
      list-style: none;
      padding-left: 20px;

      li.tag-item {
        margin-right: 10px;
        line-height: normal;

        span.placeholder {
          background-color: #ffffff;
          font-size: 0em;
          &:after {
            font-size: 0.9rem;
            font-weight: normal;
            border: 1px solid #888;
            border-radius: 3px;
            background-color: #fff;
            padding: 3px 8px;
          }
          &:hover,
          &:focus {
            &:after {
              background-color: #eee;
              border: 1px solid #666;
            }
          }
        }
        span.tag {
          display: none;
        }
      }
    }
  }

  .dropdown-content {
    display: flex;
    flex-direction: column;
    align-items: stretch;

    padding: 5px;
    padding-top: 8px;

    input.search {
      padding: 3px 8px;
      margin-left: 15px;
      margin-right: 5px;
      border: 1px solid #888;
      border-radius: 3px;
      font-size: 1em;
      &:focus {
        outline: none;
      }
    }

    ul.root {
      margin: 5px 0px;
      padding-left: 15px;
      font-weight: normal;
      font-size: 1em;

      li.node {
        i.toggle {
          font-family: monospace;
          font-size: 1.25em;
          font-style: normal;
          font-weight: 500;
          &:hover {
            color: #888888;
          }

          white-space: pre;
          margin-right: 4px;
          outline: none;

          cursor:pointer &:after {
            content: ' ';
          }
          &.collapsed:after {
            content: '+';
          }
          &.expanded:after {
            content: '-';
          }
        }

        label {
          .checkbox-item,
          .radio-item {
            vertical-align: middle;
            margin: 0 4px 0 0;
            &.simple-select {
              display: none;
            }
          }
        }
      }
    }
  }
`;

export const CoordForm = styled.form`
  display: flex;
  flex-direction: row;
  align-items: center;

  margin-top: 3px;
  margin-left: 20px;
  padding-right: 3px;

  span {
    font-weight: normal;
    margin-right: 5px;
  }
  span.coord-prefix {
    flex-shrink: 0;
  }
  input {
    flex-grow: 1;
    max-width: 15rem;
    min-width: 3rem;
  }
  button {
    flex-grow: 1;
  }
`;

export const ValidationInput = styled.input`
  border: 1px solid ${({ invalid }) => (invalid ? '#dc3545' : '#aaa')};
  &:focus {
    border: 2px solid ${({ invalid }) => (invalid ? '#dc3545' : '#aaa')};
  }
`;
ValidationInput.defaultProps = {
  invalid: false,
};

export const InvalidText = styled.span`
  margin-left: 20px;
  font-size: 0.9em;
  font-weight: normal;
  line-height: normal;
  color: #dc3545;
`;

export const RangesText = styled.span`
  margin-left: 20px;
  font-size: 0.9em;
  font-weight: normal;
  line-height: normal;
  color: #888;
`;

export const DomainSelect = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
`;
