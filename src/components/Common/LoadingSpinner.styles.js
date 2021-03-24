import styled, { keyframes } from 'styled-components';

// export const smallBoxAnim = keyframes`
//     0% {transform: scale(0.2);}
//     100% {transform: scale(0.75);}
// `;

const loaderAnim = keyframes`
    0% {transform: rotate(0deg);}
    100% {transform: rotate(90deg);}
`;

export const HollowLoader = styled.div`
  width: ${({ size }) => size};
  height: ${({ size }) => size};
  animation: ${loaderAnim} ${({ timing }) => timing} infinite ease-in-out;
  outline: 1px solid transparent;
  visibility: ${({ visible }) => (visible ? 'unset' : 'hidden')};
`;

export const LargeBox = styled.div`
  height: ${({ size }) => size};
  width: ${({ size }) => size};
  background-color: ${({ color }) => color};
  outline: 1px solid transparent;
  position: fixed;
`;

// export const SmallBox = styled.div`
//   height: ${hollowBoxSize};
//   width: ${hollowBoxSize};
//   background-color: ${hollowDark};
//   position: fixed;
//   z-index: 1;
//   outline: 1px solid transparent;
//   animation: ${smallBoxAnim} ${hollowTiming} alternate infinite ease-in-out;
// `;
