import React, { useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
// import { useStores } from '../../stores/connect';

import { PLOT_DOWNLOAD_OPTIONS } from '../../constants/defs.json';

import VegaEmbed from '../../react_vega/VegaEmbed';
import { PlotOptions } from './Plot.styles';
import DropdownButton from '../Buttons/DropdownButton';
import ExternalLink from '../Common/ExternalLink';

import initialSpec from '../../vega_specs/global_seq_v2.vg.json';

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

const DOWNLOAD_OPTIONS = {
  DOWNLOAD_SEQUENCE_DATA: 'Download global sequencing data',
  DOWNLOAD_CASE_DATA: 'Download global confirmed cases',
  DOWNLOAD_TURNAROUND: 'Download sequence submission turnaround data',
};

const GlobalSeqPlot = ({ width }) => {
  const vegaRef = useRef();
  const hiddenLink = useRef();
  // const { dataStore } = useStores();

  const handleDownloadSelect = (option) => {
    if (option === DOWNLOAD_OPTIONS.DOWNLOAD_SEQUENCE_DATA) {
      hiddenLink.current.href =
        'https://storage.googleapis.com/ve-public/new_global_data/sequences_per_month.json';
      hiddenLink.current.click();
    } else if (option === DOWNLOAD_OPTIONS.DOWNLOAD_CASE_DATA) {
      hiddenLink.current.href =
        'https://storage.googleapis.com/ve-public/new_global_data/case_count.json';
      hiddenLink.current.click();
    } else if (option === DOWNLOAD_OPTIONS.DOWNLOAD_TURNAROUND) {
      hiddenLink.current.href =
        'https://storage.googleapis.com/ve-public/new_global_data/turnaround_per_month.json';
      hiddenLink.current.click();
    } else if (option === PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_PNG) {
      vegaRef.current.downloadImage('png', 'vega-export.png', 1);
    } else if (option === PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_PNG_2X) {
      vegaRef.current.downloadImage('png', 'vega-export.png', 2);
    } else if (option === PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_PNG_4X) {
      vegaRef.current.downloadImage('png', 'vega-export.png', 4);
    } else if (option === PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_SVG) {
      vegaRef.current.downloadImage('svg', 'vega-export.svg');
    }
  };

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
        <DropdownButton
          text={'Download'}
          options={[
            DOWNLOAD_OPTIONS.DOWNLOAD_SEQUENCE_DATA,
            DOWNLOAD_OPTIONS.DOWNLOAD_CASE_DATA,
            DOWNLOAD_OPTIONS.DOWNLOAD_TURNAROUND,
            PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_PNG,
            PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_PNG_2X,
            PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_PNG_4X,
            PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_SVG,
          ]}
          onSelect={handleDownloadSelect}
        />
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
GlobalSeqPlot.propTypes = {
  width: PropTypes.number,
};
GlobalSeqPlot.defaultProps = {
  width: 1000,
};

export default GlobalSeqPlot;
