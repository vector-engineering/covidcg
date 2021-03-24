import React from 'react';

import {
  NothingSelectedContainer,
  TheText,
  Title,
} from './NothingSelected.styles';

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
