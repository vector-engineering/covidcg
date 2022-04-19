import React, { useRef, useState } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';

import VegaEmbed from '../../react_vega/VegaEmbed';

import initialSpec from '../../vega_specs/spinning_globe_seq.vg.json';

const SpinningGlobeSequences = () => {

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