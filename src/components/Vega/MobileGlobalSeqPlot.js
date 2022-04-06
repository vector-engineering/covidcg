import React, { useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
// import { useStores } from '../../stores/connect';

import VegaEmbed from '../../react_vega/VegaEmbed';
import { PlotOptions } from './Plot.styles';
import DropdownButton from '../Buttons/DropdownButton';
import ExternalLink from '../Common/ExternalLink';

import initialSpec from '../../vega_specs/mobile_global_seq_map.vg.json';

const PlotContainer = styled.div``;

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

const HelpText = styled.p`
  margin: 0px 20px;
  margin-bottom: 5px;

  font-weight: normal;
  font-size: 14px;
  line-height: normal;
`;

const MobileGlobalSeqPlot = ({ width }) => {
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
