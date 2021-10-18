import React, { useEffect } from 'react';
import styled from 'styled-components';
import ReactTooltip from 'react-tooltip';
import PropTypes from 'prop-types';

const QuestionButtonContainer = ({ ...props }) => {
  useEffect(() => {
    // This prevents every single QuestionButton from rebuilding if they aren't visible.
    if (props.rebuildAfterMount) {
      ReactTooltip.rebuild();
    }
  });
  return <span {...props}>?</span>;
};

const QuestionButton = styled(QuestionButtonContainer)`
  font-size: 0.9em;
  line-height: normal;
  font-weight: normal;

  margin-left: 8px;
  padding: 2px 5px;
  background-color: #fff;
  border: 1px solid #ccc;
  border-radius: 2px;
  &:hover {
    background-color: #f8f8f8;
  }
`;

QuestionButtonContainer.defaultProps = {
  rebuild: ''
};

QuestionButtonContainer.propTypes = {
  rebuild: PropTypes.bool
};

export default QuestionButton;