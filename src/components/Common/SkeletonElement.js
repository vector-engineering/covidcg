import React, { memo } from 'react';
import PropTypes from 'prop-types';

import { SkeletonContainer } from './SkeletonElement.styles';

// eslint-disable-next-line react/prop-types
const SkeletonElement = memo(({ width, height, delay, style, children }) => {
  // console.log('render: ', delay);
  return (
    <SkeletonContainer
      width={width}
      height={height}
      delay={delay}
      style={style}
    >
      {children}
    </SkeletonContainer>
  );
});
SkeletonElement.propTypes = {
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node,
  ]),
  width: PropTypes.string,
  height: PropTypes.number.isRequired,
  delay: PropTypes.number.isRequired,
  style: PropTypes.object,
};
SkeletonElement.defaultProps = {
  width: '100%',
  style: {},
};

SkeletonElement.displayName = 'SkeletonElement';

export default SkeletonElement;
