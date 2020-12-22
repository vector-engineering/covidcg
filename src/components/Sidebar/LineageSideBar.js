import React from 'react';
import styled from 'styled-components';
import VegaLegend from '../Legend/VegaLegend';

const Container = styled.div`
  width: 90px;
  height: 100vh;
  overflow-y: scroll;
  border-left: 1px #eaeaea solid;
`;

const LineageSidebar = () => {
  return (
    <Container>
      <VegaLegend />
    </Container>
  );
};

export default LineageSidebar;
