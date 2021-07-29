import styled from 'styled-components';

import DropdownTreeSelect from 'react-dropdown-tree-select';

export const ContainerDiv = styled.div`
  width: 100%;
  height: 100%;

  margin-top: 2px;
  padding-top: 8px;

  display: flex;
  flex-direction: column;
  // overflow-y: hidden;
  overflow-y: scroll;
`;

export const DropdownHeader = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;

  margin-left: 15px;
`;

export const Title = styled.span`
  font-size: 1rem;
  font-weight: 500;
`;

export const UnselectButton = styled.button`
  display: ${({ show }) => (show ? 'block' : 'none')};
  margin-left: 10px;
`;
UnselectButton.defaultProps = {
  show: true,
};

export const SelectedLocationsContainer = styled.div`
  display: flex;
  flex-direction: row;
  flex-wrap: wrap;
  padding: 5px 0px;
`;

export const StyledDropdownTreeSelect = styled(DropdownTreeSelect)`
  margin-top: 3px;
  overflow-y: hidden;

  span.placeholder {
    background-color: #ffffff;
    font-size: 0em;
  }

  .node > label {
    cursor: pointer;
    margin-left: 2px;
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
    height: 100%;
    flex-direction: column;
    display: flex;
    // overflow: hidden;
    align-items: stretch;

    a.dropdown-trigger {
      display: none;
    }

    .dropdown-content {
      position: relative;
      height: calc(100vh - 220px);
      padding: 4px;
      padding-left: 15px;
      padding-right: 15px;
      background-color: #f8f8f8;
      overflow-y: scroll;

      input.search {
        font-size: 1em;
        padding: 5px 8px;
        width: calc(100% - 20px);
        border: 1px solid #aaa;
        border-radius: 3px;
        background-color: #ffffff;
        outline: none;
      }

      ul.root {
        margin-top: 5px;
        padding: 0;
        flex-direction: column;
        display: flex;
        overflow-y: scroll;

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
        .infinite-scroll-component {
          overflow-x: hidden;
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

  .select-all-children {
    font-style: normal;
    font-weight: normal;

    border: 1px solid #ccc;
    padding: 1px 4px;
    border-radius: 3px;
    background-color: #fff;

    margin-left: 5px;

    cursor: pointer;

    &:after {
      content: '↳';
    }

    &:hover,
    &:focus {
      border-color: #888;
      background-color: #f8f8f8;
    }

    &:active {
      border-color: #000;
      background-color: #eee;
    }
  }
`;
