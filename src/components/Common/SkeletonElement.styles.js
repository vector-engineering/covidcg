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

export const SkeletonContainer = styled.div`
  height: ${({ height }) => height}px;
  width: ${({ width }) => width};
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
