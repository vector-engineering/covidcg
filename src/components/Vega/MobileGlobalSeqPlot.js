import React, { useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
// import { useStores } from '../../stores/connect';

import VegaEmbed from '../../react_vega/VegaEmbed';
import { PlotOptions } from './Plot.styles';
import DropdownButton from '../Buttons/DropdownButton';
import ExternalLink from '../Common/ExternalLink';

import initialSpec from '../../vega_specs/mobile_global_seq_map.vg.json';

const MobileGlobalSeqPlot = () => {
    const vegaRef = useRef();
    // const { dataStore } = useStores();

  return (
            <VegaEmbed
                ref={vegaRef}
                spec={initialSpec}
                width={window.width}
                actions={false}
            />
    );
};
MobileGlobalSeqPlot.propTypes = {
    width: PropTypes.number,
};
MobileGlobalSeqPlot.defaultProps = {
    width: 375,
};

export default MobileGlobalSeqPlot;