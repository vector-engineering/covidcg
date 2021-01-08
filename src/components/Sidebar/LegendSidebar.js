import React from 'react';
import styled from 'styled-components';
import LegendContainer from '../Legend/LegendContainer';

const Container = styled.div`
  width: 150px;
  height: 100vh;
  overflow-y: hidden;
  overflow-x: hidden;
  border-left: 1px #eaeaea solid;
`;

const LegendSidebar = () => {
  return (
    <Container>
      <LegendContainer />
    </Container>
  );
};

export default LegendSidebar;
