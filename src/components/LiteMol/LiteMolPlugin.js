import React, { useEffect } from 'react';
import styled from 'styled-components';
import PropTypes from 'prop-types';

import LiteMol from 'litemol';
const Bootstrap = LiteMol.Bootstrap;

const PluginContainer = styled.div`
  height: ${({ height }) => height}px;
  .lm-plugin {
    height: ${({ height }) => height}px;
  }
`;

const LiteMolPlugin = React.forwardRef(({ setPlugin, height }, ref) => {
  useEffect(() => {
    const _plugin = LiteMol.Plugin.create({
      target: ref.current,
      layoutState: {
        hideControls: true,
        isExpanded: false,
        collapsedControlsLayout:
          Bootstrap.Components.CollapsedControlsLayout.Landscape,
      },
      viewportBackground: '#AAAAAA',
    });
    setPlugin(_plugin);
  }, []);

  return <PluginContainer height={height} ref={ref} />;
});
LiteMolPlugin.propTypes = {
  setPlugin: PropTypes.func,
  height: PropTypes.number,
};
LiteMolPlugin.defaultProps = {
  height: 500,
};

export default LiteMolPlugin;
