import React, { useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';

import VegaEmbed from '../../react_vega/VegaEmbed';
import initialSpec from '../../vega_specs/map_combined_standalone.vg.json';

const PlotContainer = styled.div``;

const SequencingMapPlot = ({ width }) => {
  const vegaRef = useRef();

  return (
    <PlotContainer>
      <VegaEmbed
        ref={vegaRef}
        spec={initialSpec}
        width={width}
        actions={false}
      />
    </PlotContainer>
  );
};
SequencingMapPlot.propTypes = {
  width: PropTypes.number,
};
SequencingMapPlot.defaultProps = {
  width: 100,
};

export default SequencingMapPlot;
