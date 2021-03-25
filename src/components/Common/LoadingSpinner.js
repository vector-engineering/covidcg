import React from 'react';
import PropTypes from 'prop-types';
import styled, { keyframes } from 'styled-components';

// const smallBoxAnim = keyframes`
//     0% {transform: scale(0.2);}
//     100% {transform: scale(0.75);}
// `;

const loaderAnim = keyframes`
    0% {transform: rotate(0deg);}
    100% {transform: rotate(90deg);}
`;

const HollowLoader = styled.div`
  width: ${({ size }) => size};
  height: ${({ size }) => size};
  animation: ${loaderAnim} ${({ timing }) => timing} infinite ease-in-out;
  outline: 1px solid transparent;
  visibility: ${({ visible }) => (visible ? 'unset' : 'hidden')};
`;

const LargeBox = styled.div`
  height: ${({ size }) => size};
  width: ${({ size }) => size};
  background-color: ${({ color }) => color};
  outline: 1px solid transparent;
  position: fixed;
`;

// const SmallBox = styled.div`
//   height: ${hollowBoxSize};
//   width: ${hollowBoxSize};
//   background-color: ${hollowDark};
//   position: fixed;
//   z-index: 1;
//   outline: 1px solid transparent;
//   animation: ${smallBoxAnim} ${hollowTiming} alternate infinite ease-in-out;
// `;

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
