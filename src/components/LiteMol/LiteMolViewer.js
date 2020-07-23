import React, { useLayoutEffect, useRef } from 'react';
import styled from 'styled-components';
import LiteMol from 'litemol';
import './../../styles/litemol.min.css';

const Container = styled.div`
  padding-top: 20px;
  width: 100%;
  height: 100%;
`;

const LiteMolViewer = React.memo(() => {
  const target = useRef();
  const plugin = useRef();

  useLayoutEffect(() => {
    plugin.current = LiteMol.Plugin.create({
      target: target.current,
      layoutState: {
        hideControls: true,
        collapsedControlsLayout:
          LiteMol.Bootstrap.Components.CollapsedControlsLayout.Landscape,
      },
      viewportBackground: '#F1F1F1',
    });
    plugin.current.loadMolecule({
      id: '1tqn',
      url: 'https://files.rcsb.org/download/6X2A.pdb',
      format: 'pdb', // default
    });
  });

  return (
    <Container>
      <div
        style={{ width: '100%', height: '100%' }}
        id={`${Math.random()}-litemolViewer`}
        ref={target}
      />
    </Container>
  );
});

export default LiteMolViewer;
