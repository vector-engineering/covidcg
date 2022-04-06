import React, { useRef, useState } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';

import VegaEmbed from '../../react_vega/VegaEmbed';

import initialSpec from '../../vega_specs/spinning_globe_seq.vg.json';

import scores from '../../vega_specs/map_scores';
const PlotContainer = styled.div`
width: 100%;`;

const PlotTitle = styled.div`
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  justify-content: center;

  margin-right: 10px;
  padding-right: 10px;
  padding-left: 18px;

  line-height: normal;

  .title {
    font-size: 1.25rem;
  }
  .subtitle {
    font-size: 0.9em;
    font-weight: normal;
  }
`;

// TODO: in-spec: stop scrolling when interacted with.

const SpinningGlobeSequences = ({ width }) => {
    const {
        dataStore
    } = useStores();

    // Generates data to be passed into the globe.
    const genData = () => {
        let d = dataStore;
        //console.log(JSON.stringify(dataStore));

    };

    const vegaRef = useRef();
    // const { dataStore } = useStores();

    genData();
    return (
            <VegaEmbed
                ref={vegaRef}
                spec={initialSpec}
                width={window.width}
                actions={false}
            />
    );
};
SpinningGlobeSequences.propTypes = {
    width: PropTypes.number,
};
SpinningGlobeSequences.defaultProps = {
    width: 375,
};

export default SpinningGlobeSequences;
