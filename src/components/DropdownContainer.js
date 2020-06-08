import React, { Component } from 'react';
import PropTypes from 'prop-types';
import DropdownTreeSelect from 'react-dropdown-tree-select';
import styled from 'styled-components';
import _ from 'underscore';

const StyledDropdownTreeSelect = styled(DropdownTreeSelect)`
  .tag {
    background-color: #f4f4f4;
    border: 1px solid #e9e9e9;
    padding: 2px 0 2px 2px;
    border-radius: 2px;
    display: inline-block;
    &:focus-within {
      background-color: #e9e9e9;
      border-color: #a0a0a0;
    }
  }
  .tag-remove {
    color: #a0a0a0;
    font-size: 75%;
    line-height: 100%;
    cursor: pointer;
    background-color: transparent;
    border: none;
    outline: none;
    &.disabled,
    &.readOnly {
      cursor: not-allowed;
    }
    &:focus {
      color: #3c3c3c;
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
      width: 100%;
      background-color: #ffffff;
      border: none;
      border-top: 1px solid #aaa;
      border-bottom: 1px solid #aaa;
      box-shadow: 0px 0px 3px #aaa;
      padding: 8px 0px;
      line-height: 20px;
      max-height: 200px;
      display: inline-block;
      overflow: auto;

      &:focus {
        //outline: none;
      }
      &.arrow {
        cursor: pointer;
        &.bottom:after {
          content: '\\25BC';
          vertical-align: middle;
          color: #3c3c3c;
          margin-right: 2px;
        }
        &.top:after {
          content: '\\25B2';
          vertical-align: middle;
          color: #3c3c3c;
          margin-right: 2px;
        }
        &.disabled {
          cursor: not-allowed;
        }
        &.disabled.bottom:after {
          color: #b9b9b9;
        }
      }
      . ul {
        margin-left: 5px;

        li.tag-item {
          span {
            padding: 2px 5px;
            background-color: #ffffff;
            border: 1px solid #aaaaaa;
            box-shadow: 0px 0px 1px 2px #eee;
            border-radius: 5px;

            button.tag-remove {
              font-size: 1em;
              color: #e67272;
              font-weight: 800;
            }

            &.placeholder {
              padding: 5px 5px;
            }
          }
        }
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
