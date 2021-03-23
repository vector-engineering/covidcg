import styled from 'styled-components';

export const NothingSelectedContainer = styled.div`
  width: 100%;
  height: 100%;
  position: absolute;
  background-color: white;
  top: 0;
  left: 0;
  z-index: 9999;
  display: flex;
  align-items: center;
  justify-content: flex-start;
  flex-direction: column;
  padding-top: 20%;
`;

export const TheText = styled.div`
  font-weight: 300;
  max-width: 400px;
  margin-bottom: 40px;
`;

export const Title = styled.div`
  font-weight: 500;
  font-size: 32px;
  line-height: 32px;
`;
