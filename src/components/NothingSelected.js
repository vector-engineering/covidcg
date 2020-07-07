import React from 'react';
import styled from 'styled-components';

const NothingSelectedContainer = styled.div`
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

const TheText = styled.div`
  font-weight: 300;
  max-width: 400px;
  margin-bottom: 40px;
`;

const Title = styled.div`
  font-weight: 500;
  font-size: 32px;
  line-height: 32px;
`;

const NothingSelected = () => {
  return (
    <NothingSelectedContainer>
      <TheText>
        <Title>Compare lineages by location in seconds</Title>
        <br />
        No sequences selected, select a location on the sidebar to begin
      </TheText>
      <img
        style={{ width: '400px', height: 'auto' }}
        src="https://i.imgur.com/3plI4KX.png"
      />
    </NothingSelectedContainer>
  );
};

export default NothingSelected;
