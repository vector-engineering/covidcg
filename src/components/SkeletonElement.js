import React, { memo } from 'react';
import PropTypes from 'prop-types';
import styled, { keyframes } from 'styled-components';

const loadingAnim = keyframes`
    0%,
    15%,
    85%,
    100% {
      opacity: 1;
    }
    60% {
      opacity: 0.2;
    }
  
`;

const SkeletonContainer = styled.div`
  height: ${({ height }) => height};
  width: 100%;
  background-color: #eee;
  animation-name: ${loadingAnim};
  animation-duration: 2250ms;
  transition-timing-function: ease-in-out;
  animation-iteration-count: infinite;
  animation-delay: ${({ delay }) => `${delay * 95}ms`};
  display: flex;
  justify-content: center;
  align-items: center;
`;

const SkeletonElement = memo(({ height, delay, children }) => {
  // console.log('render: ', delay);
  return (
    <SkeletonContainer height={height} delay={delay}>
      {children}
    </SkeletonContainer>
  );
});
SkeletonElement.propTypes = {
  height: PropTypes.number.isRequired,
  delay: PropTypes.number.isRequired,
};

SkeletonElement.displayName = 'SkeletonElement';

export default SkeletonElement;
