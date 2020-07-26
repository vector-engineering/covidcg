import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import VegaEmbed from '../../react_vega/VegaEmbed';
import initialSpec from '../../vega/location_date.vg.json';

const PlotContainer = styled.div``;

const LocationDatePlot = observer(({ width }) => {
  const vegaRef = useRef();
  const { covidStore } = useStores();

  const [state, setState] = useState({
    data: {
      location_data: [],
      selectedGroups: [],
    },
    spec: JSON.parse(JSON.stringify(initialSpec)),
    hoverGroup: {},
  });

  useEffect(() => {
    setState({
      ...state,
      data: {
        location_data: JSON.parse(JSON.stringify(covidStore.aggLocationData)),
        selectedGroups: JSON.parse(JSON.stringify(covidStore.selectedGroups)),
      },
    });
  }, [covidStore.aggLocationData, covidStore.selectedGroups]);

  // useEffect(() => {
  //   let spec = JSON.parse(JSON.stringify(state.spec));

  //   // Set the width manually
  //   // spec['width'] = width;
  //   // spec['height'] = 300;

  //   setState({ ...state, spec });
  // });

  return (
    <PlotContainer>
      <div style={{ width: `${width}` }}>
        <VegaEmbed
          ref={vegaRef}
          data={state.data}
          spec={state.spec}
          width={width}
          actions={false}
        />
      </div>
    </PlotContainer>
  );
});
LocationDatePlot.propTypes = {
  width: PropTypes.number,
};
LocationDatePlot.defaultProps = {
  width: 100,
};

export default LocationDatePlot;
