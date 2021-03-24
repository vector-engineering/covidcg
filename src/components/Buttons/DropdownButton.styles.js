import styled from 'styled-components';

export const DropdownContainer = styled.div`
  position: relative;

  button {
    padding-right: 18px;
  }

  .caret {
    position: relative;
    margin-left: 5px;
  }

  .caret:after {
    content: '';
    position: absolute;
    left: 1px;
    top: 5px;
    border-top: 5px solid #eeeeee;
    border-left: 5px solid transparent;
    border-right: 5px solid transparent;
  }
`;

export const DropdownMenu = styled.div`
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  justify-content: center;

  position: absolute;
  top: 28px;
  ${({ direction }) => (direction === 'left' ? 'left: 0' : 'right: 0')};
  z-index: 10;

  min-width: 10rem;
  padding: 5px 0px;

  background-color: #fff;
  border: 1px solid #aaa;
  border-radius: 3px;

  overflow: hidden;
`;

DropdownMenu.defaultProps = {
  direction: 'right',
};

export const DropdownItem = styled.a`
  display: block;
  width: 100%;
  padding: 0.25rem 1.5rem;
  clear: both;
  font-weight: 400;
  color: #212529;
  text-align: inherit;
  white-space: nowrap;
  background-color: #fff;
  border: 0;
  text-decoration: none;

  &:hover {
    background-color: #f8f9fa;
  }
  &:active {
    color: #ffffff;
    background-color: #34d058;
  }
`;
