import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';

const AckCellContainer = styled.div`
  white-space: pre-line;
  line-height: initial;
  padding: 4px 0px;
`;

const AckCell = ({ text }) => {
  return <AckCellContainer>{text}</AckCellContainer>;
};
AckCell.displayName = 'AckCell';
AckCell.propTypes = {
  text: PropTypes.string,
};
AckCell.defaultProps = {
  text: '',
};

export default AckCell;
