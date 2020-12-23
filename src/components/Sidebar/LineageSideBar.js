import React from 'react';
import styled from 'styled-components';
import LegendContainer from '../Legend/LegendContainer';

const Container = styled.div`
  width: 90px;
  height: 100vh;
  overflow-y: scroll;
  border-left: 1px #eaeaea solid;
`;

const LineageSidebar = () => {
  return (
    <Container>
      <LegendContainer />
    </Container>
  );
};

export default LineageSidebar;
