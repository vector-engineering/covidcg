import React, { Component } from 'react';
import PropTypes from 'prop-types';
import DropdownTreeSelect from 'react-dropdown-tree-select';
import styled from 'styled-components';
import _ from 'underscore';

const StyledDropdownTreeSelect = styled(DropdownTreeSelect)`
  ul.tag-list {
    margin-left: 5px;
    li:first-child {
      span.placeholder:after {
        content: 'Viewing all sequences';
        font-size: 0.9rem;
        font-weight: normal;
      }
    }
  }
  .tag {
    background-color: #f4f4f4;
    border: 1px solid #e9e9e9;
    padding: 3px 6px;
    border-radius: 2px;
    display: inline-block;
    font-weight: normal;
    &:focus-within {
      background-color: #e9e9e9;
      border-color: #a0a0a0;
    }
  }
  span.placeholder {
    background-color: #ffffff;
    font-size: 0em;
  }
  .tag-remove {
    color: #a0a0a0;
    font-size: 1.25em;
    line-height: 100%;
    cursor: pointer;
    background-color: transparent;
    border: none;
    outline: none;
    &:hover,
    &:focus {
      color: #ff5555;
    }
    &.disabled,
    &.readOnly {
      cursor: not-allowed;
    }
  }

  .node > label {
    cursor: pointer;
    margin-left: 2px;
  }

  .tag-list {
    display: inline;
    padding: 0;
    margin: 0;
  }

  .tag-item {
    display: inline-block;
    margin: 4px;

    .search {
      border: none;
      border-bottom: 1px solid #ccc;
      outline: none;
    }
    &:last-child {
      margin-right: 4px;
    }
  }

  .node {
    list-style: none;
    white-space: nowrap;
    padding: 4px;
    &.leaf.collapsed {
      display: none;
    }
    &.disabled > * {
      color: gray;
      cursor: not-allowed;
    }
    &.match-in-children.hide {
      .node-label {
        opacity: 0.5;
      }
    }
    &.focused {
      background-color: #f4f4f4;
    }
    .toggle {
      white-space: pre;
      margin-right: 4px;
      cursor: pointer;
    }
  }
  .toggle {
    white-space: pre;
    margin-right: 4px;
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

  .searchModeOn {
    .toggle {
      display: none;
    }
  }

  .checkbox-item,
  .radio-item {
    vertical-align: middle;
    margin: 0 4px 0 0;
    &.simple-select {
      display: none;
    }
  }

  .hide:not(.match-in-children) {
    display: none;
  }

  .dropdown {
    width: 100%;
    display: block;
    position: relative;

    a.dropdown-trigger {
      width: calc(100% - 16px);
      background-color: #ffffff;
      border: none;
      border-top: 1px solid #aaa;
      border-bottom: 1px solid #aaa;
      box-shadow: 0px 0px 2px #aaa;
      padding: 8px;
      line-height: 20px;
      max-height: 200px;
      display: inline-block;
      overflow: auto;

      &:focus {
        outline: none;
      }
    }

    .dropdown-content {
      position: relative;
      //width: calc(100% - 10px);
      padding: 4px;
      padding-left: 8px;
      background-color: #f8f8f8;
      z-index: 1;

      input.search {
        font-size: 1em;
        padding: 5px 8px;
        width: calc(100% - 20px);
        border: 1px solid #aaa;
        border-radius: 5px;
        background-color: #ffffff;
        outline: none;
      }

      ul.root {
        margin-top: 5px;
        padding: 0;

        i.toggle {
          font-family: monospace;
          font-size: 1.25em;
          font-style: normal;
          font-weight: 800;
          &:hover {
            color: #888888;
          }
          &:focus {
            outline: none;
          }
        }
      }
    }
  }

  // Custom node styling
  .fa.fa-info {
    margin-left: 0.5em;
    font-weight: normal;
    font-style: normal;
  }
`;

class DropdownContainer extends Component {
  // I forget why this component needs a wrapper with it's own state management
  // but for some reason it does. If I remove this stuff the selections don't work
  // anymore.
  constructor(props) {
    super(props);
    this.state = {
      data: props.data,
    };
  }

  UNSAFE_componentWillReceiveProps = (nextProps) => {
    if (!_.isEqual(nextProps.data, this.state.data)) {
      this.setState({ data: nextProps.data });
    }
  };

  shouldComponentUpdate = (nextProps) => {
    return !_.isEqual(nextProps.data, this.state.data);
  };

  render() {
    const { data, ...rest } = this.props;
    data;
    return <StyledDropdownTreeSelect data={this.state.data} {...rest} />;
  }
}

DropdownContainer.propTypes = {
  data: PropTypes.oneOfType([PropTypes.object, PropTypes.array]).isRequired,
};

export default DropdownContainer;
