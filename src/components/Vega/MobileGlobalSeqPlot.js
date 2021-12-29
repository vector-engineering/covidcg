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
        <PlotContainer>
            <PlotOptions>
                <PlotTitle>
                    <span className="title">Global Sequencing Effort</span>
                </PlotTitle>
                <div className="spacer"></div>
            </PlotOptions>
            <HelpText>
                The number of genomic sequence and associate data are shared via the
                GISAID Initiative (
                <ExternalLink href="https://doi.org/10.1002/gch2.1018">
                    Elbe et al, 2017, <i>Wiley Global Challenges</i>
                </ExternalLink>
                ) and case data is obtained from{' '}
                <ExternalLink href="https://github.com/CSSEGISandData/COVID-19">
                    JHU CSSE COVID-19 Data
                </ExternalLink>{' '}
                (
                <ExternalLink href="https://doi.org/10.1016/S1473-3099(20)30120-1">
                    Dong et al, 2020, <i>Lancet Inf Dis.</i>
                </ExternalLink>
                ). Data from cases and sequences are grouped by month to reduce noise.
                Regions with &gt;20 sequences per 1000 cases are colored the same in the
                left map.
            </HelpText>
            <HelpText>
                Hover over points in the top-left plot, or over countries on the map, to
                highlight them and display them on the bottom row of plots. Click on a
                point, or double-click a country, to track that selection and
                persistently display them on the bottom row of plots.
            </HelpText>
            <HelpText>
                Click and drag to move the map, and scroll up or down to zoom in or out.
            </HelpText>
            <HelpText>
                <i>
                    *Please note that we are currently working on cleaning the metadata
                    delineating territories and countries
                </i>
            </HelpText>
            <VegaEmbed
                ref={vegaRef}
                spec={initialSpec}
                width={width}
                actions={false}
            />
        </PlotContainer>
    );
};
MobileGlobalSeqPlot.propTypes = {
    width: PropTypes.number,
};
MobileGlobalSeqPlot.defaultProps = {
    width: 400,
};

export default MobileGlobalSeqPlot;
