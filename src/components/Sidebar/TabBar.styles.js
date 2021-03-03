import styled from 'styled-components';

export const TabBarContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: stretch;
  flex-shrink: 0;

  width: 100%;
  min-height: 30px;
  height=${({ height }) => height}px;
  background-color: #eee;
`;

export const TabBarList = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;
  justify-content: flex-start;
  flex-grow: 1;
  border-bottom: 1px solid #ccc;
`;

export const TabItem = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  border-bottom: 1px solid #ccc;

  height: 30px;

  a.tab-link {
    display: flex;
    flex-direction: column;
    align-items: stretch;

    flex-grow: 1;
    text-decoration: none;
    color: ${({ active }) => (active ? '#000' : '#888')};
    background-color: ${({ active }) => (active ? '#fff' : 'transparent')};
    transition: 0.1s all ease-in-out;

    padding: 4px 10px;

    &:hover {
      color: #666;
      background-color: ${({ active }) => (active ? '#fff' : '#f8f8f8')};
    }
    &:active {
      color: #888;
    }

    span {
      // border-right: ${({ active }) => (active ? 'none' : '1px solid #ccc')};
    }
  }
`;
TabItem.defaultProps = {
  active: false,
};

export const DropdownTab = styled.button`
  display: flex;
  flex-direction: row;
  align-items: center;

  width: 100%;
  height: 30px;
  padding: 0px 15px;
  padding-right: 20px;

  border: none;
  outline: none;

  background-color: ${({ active }) => (active ? '#fff' : 'transparent')};
  line-height: normal;
  text-align: center;
  font: 14px 'Helvetica Neue', Helvetica, Arial, sans-serif;
  font-weight: 500;
  color: ${({ active }) => (active ? '#000' : '#888')};

  transition: 0.1s all ease-in-out;

  .caret {
    height: 5px;

    &:after {
      left: 1px;
      top: 0px;
      border-top: 5px solid #888;
    }
  }
`;
DropdownTab.defaultProps = {
  active: false,
};
