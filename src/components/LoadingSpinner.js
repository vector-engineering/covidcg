import React from 'react';
import styled, { keyframes } from 'styled-components';

const hollowBoxSize = '3em';
const hollowLight = '#34495e';
// const hollowDark = '#34495e';
const hollowTiming = '1125ms';

// const smallBoxAnim = keyframes`
//     0% {transform: scale(0.2);}
//     100% {transform: scale(0.75);}
// `;

const loaderAnim = keyframes`
    0% {transform: rotate(0deg);}
    100% {transform: rotate(90deg);}
`;

const HollowLoader = styled.div`
  width: ${hollowBoxSize};
  height: ${hollowBoxSize};
  animation: ${loaderAnim} ${hollowTiming} infinite ease-in-out;
  outline: 1px solid transparent;
`;

const LargeBox = styled.div`
  height: ${hollowBoxSize};
  width: ${hollowBoxSize};
  background-color: ${hollowLight};
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

const LoadingSpinner = () => {
  return (
    <HollowLoader>
      <LargeBox />
    </HollowLoader>
  );
};

export default LoadingSpinner;
