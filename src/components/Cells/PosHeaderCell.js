import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';

const PosHeaderCellContainer = styled.div`
  transform: rotate(-55deg) translate(8px, 0px);
`;

const PosHeaderCell = ({ pos }) => {
  return <PosHeaderCellContainer>{pos}</PosHeaderCellContainer>;
};

PosHeaderCell.propTypes = {
  pos: PropTypes.oneOfType([PropTypes.string, PropTypes.number]),
};
PosHeaderCell.defaultProps = {
  pos: 0,
};

export default PosHeaderCell;
