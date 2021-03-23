import React from 'react';
import PropTypes from 'prop-types';

import { HollowLoader, LargeBox } from './LoadingSpinner.styles';

const LoadingSpinner = ({ size, color, timing, visible }) => {
  return (
    <HollowLoader size={size} timing={timing} visible={visible}>
      <LargeBox size={size} color={color} />
    </HollowLoader>
  );
};
LoadingSpinner.propTypes = {
  size: PropTypes.string,
  color: PropTypes.string,
  timing: PropTypes.string,
  visible: PropTypes.bool,
};
LoadingSpinner.defaultProps = {
  size: '2em',
  color: '#ccc', // dark: #34495e
  timing: '1125ms',
  visible: true,
};

export default LoadingSpinner;
