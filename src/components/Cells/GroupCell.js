import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';

const GroupCellContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: stretch;

  margin-left: -6px;

  .stripe {
    flex-shrink: 0;
    background-color: ${({ color }) => color};
    width: 10px;
    margin-right: 5px;
  }

  .text {
    flex-grow: 1;
  }
`;
GroupCellContainer.defaultProps = {
  color: '#fff',
};

const GroupCell = ({ text, color }) => {
  return (
    <GroupCellContainer color={color}>
      <div className="stripe" />
      <span className="text">{text}</span>
    </GroupCellContainer>
  );
};
GroupCell.displayName = 'GroupCell';
GroupCell.propTypes = {
  text: PropTypes.oneOfType([PropTypes.string, PropTypes.number]),
  color: PropTypes.string,
};
GroupCell.defaultProps = {
  text: '',
  color: '#fff',
};

export default GroupCell;
