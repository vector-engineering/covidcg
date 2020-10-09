import React, { useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';

import { PLOT_DOWNLOAD_OPTIONS } from '../../constants/download';

import VegaEmbed from '../../react_vega/VegaEmbed';
import { PlotOptions } from './Plot.styles';
import DropdownButton from '../Buttons/DropdownButton';

import initialSpec from '../../vega_specs/map_combined.vg.json';
// eslint-disable-next-line import/no-unresolved
import countryScoreData from 'https://storage.googleapis.com/ve-public/country_score.json';

const PlotContainer = styled.div``;

const HelpText = styled.p`
  margin: 0px 20px;
  font-weight: normal;
  font-size: 14px;
  line-height: normal;
`;

const SequencingMapPlot = ({ width }) => {
  const vegaRef = useRef();
  const { dataStore } = useStores();

  const initialData = {
    scores: countryScoreData,
  };

  const handleDownloadSelect = (option) => {
    if (option === PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_DATA) {
      dataStore.downloadCountryScoreData();
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
      <PlotOptions>
        <HelpText>
          Click and drag to move the map. Scroll or use the mouse wheel to zoom
          in and out of the map.
        </HelpText>
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
