import React, { useRef, useState } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';

import { PLOT_DOWNLOAD_OPTIONS } from '../../constants/defs.json';

import VegaEmbed from '../../react_vega/VegaEmbed';
import { PlotOptions } from './Plot.styles';
import ExternalLink from '../Common/ExternalLink';

import initialSpec from '../../vega_specs/spinning_globe_seq.vg.json';
import { StatusText } from '../Sidebar/StatusBox.styles';

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

const HelpText = styled.p`
  margin: 0px 20px;
  margin-bottom: 5px;

  font-weight: normal;
  font-size: 14px;
  line-height: normal;
`;

const SpinningGlobeSequences = ({ width }) => {
    const {
        dataStore
    } = useStores();

    const [state, setState] = useState({
        data: {
            scores: []
        }
        }
    );

    // Generates data to be passed into the globe.
    const genData = () => {
        let d = dataStore;
        console.log(JSON.stringify(dataStore));

    };

    const vegaRef = useRef();
    const hiddenLink = useRef();
    // const { dataStore } = useStores();

    genData();
    return (
        <PlotContainer>
            <a
                ref={hiddenLink}
                href=""
                target="_blank"
                rel="noopener noreferrer"
                style={{ visibility: 'hidden' }}
            />
            <PlotOptions>
                <PlotTitle>
                    <span className="title">Global Sequencing Effort</span>
                </PlotTitle>
                <div className="spacer"></div>
            </PlotOptions>
            <VegaEmbed
                ref={vegaRef}
                spec={initialSpec}
                data={state.data}
                width={window.width}
                actions={false}
            />
        </PlotContainer>
    );
};
SpinningGlobeSequences.propTypes = {
    width: PropTypes.number,
};
SpinningGlobeSequences.defaultProps = {
    width: 1000,
};

export default SpinningGlobeSequences;
