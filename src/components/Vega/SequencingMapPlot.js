import React, { useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';

import VegaEmbed from '../../react_vega/VegaEmbed';
import initialSpec from '../../vega_specs/map_combined.vg.json';
import countryScoreData from '../../../data/country_score.json';

const PlotContainer = styled.div``;

const SequencingMapPlot = ({ width }) => {
  const vegaRef = useRef();

  const initialData = {
    scores: countryScoreData,
  };

  return (
    <PlotContainer>
      <VegaEmbed
        ref={vegaRef}
        spec={initialSpec}
        data={initialData}
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
