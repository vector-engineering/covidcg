import React, { useRef, useState } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
// import { useStores } from '../../stores/connect';

import { PLOT_DOWNLOAD_OPTIONS } from '../../constants/defs.json';

import VegaEmbed from '../../react_vega/VegaEmbed';
import { PlotOptions } from './Plot.styles';
import DropdownButton from '../Buttons/DropdownButton';
import WarningBox from '../Common/WarningBox';

import initialSpec from '../../vega_specs/surveillance.vg.json';

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

const SurveillancePlot = ({ width }) => {
  const vegaRef = useRef();
  const hiddenLink = useRef();
  // const { dataStore } = useStores();

  const [state, setState] = useState({
    showWarning: true,
  });

  const onDismissWarning = () => {
    setState({
      ...state,
      showWarning: false,
    });
  };

  const handleDownloadSelect = (option) => {
    if (option === PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_DATA) {
      hiddenLink.current.href =
        'https://storage.googleapis.com/ve-public/surveillance/group_counts.json';
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
          <span className="title">Global Lineage Surveillance</span>
        </PlotTitle>
        <div className="spacer"></div>
        <DropdownButton
          text={'Download'}
          options={[
            PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_DATA,
            PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_PNG,
            PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_PNG_2X,
            PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_PNG_4X,
            PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_SVG,
          ]}
          onSelect={handleDownloadSelect}
        />
      </PlotOptions>
      <HelpText>
        Only data from the last 90 days is shown. Please note that the most
        recent data is sparser due to lags in time between sample collection and
        submission.
      </HelpText>
      <HelpText>
        Lineages with increased prevalence in at least one region are colored
        and shown in the legend to the left. Hover over lineages in the legend,
        or near them in the plots, to highlight the lineage across all plots.
      </HelpText>
      <WarningBox show={state.showWarning} onDismiss={onDismissWarning}>
        Inconsistent sampling in the underlying data can result in missing data
        and artefacts in this visualization. Increased prevalence of lineages{' '}
        <b>does not</b>, on its own, suggest an increase in transmissibility.
        Please interpret this data with care.
      </WarningBox>
      <VegaEmbed
        ref={vegaRef}
        spec={initialSpec}
        width={width}
        actions={false}
      />
    </PlotContainer>
  );
};
SurveillancePlot.propTypes = {
  width: PropTypes.number,
};
SurveillancePlot.defaultProps = {
  width: 1000,
};

export default SurveillancePlot;
