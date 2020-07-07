import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';

const AckAuthorCellContainer = styled.div`
  white-space: pre-line;
  line-height: initial;
  padding: 4px 0px;

  a {
    text-decoration: none;
    color: #333;
  }
`;

const AckAuthorCell = ({ shortText, longText }) => {
  return (
    <AckAuthorCellContainer>
      <a title={longText}>{shortText}</a>
    </AckAuthorCellContainer>
  );
};
AckAuthorCell.displayName = 'AckAuthorCell';
AckAuthorCell.propTypes = {
  shortText: PropTypes.string.isRequired,
  longText: PropTypes.string.isRequired,
};
AckAuthorCell.defaultProps = {
  shortText: '',
  longText: '',
};

export default AckAuthorCell;
