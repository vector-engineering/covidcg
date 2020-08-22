import React, { useState, useEffect } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { toJS } from 'mobx';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import _ from 'underscore';

import VegaEmbed from '../../react_vega/VegaEmbed';

import initialSpec from '../../vega_specs/entropy.vg.json';
import { COORDINATE_MODES, DNA_OR_AA } from '../../constants/config';

const PlotContainer = styled.div``;

const EntropyPlot = observer(({ width }) => {
  const { configStore, dataStore } = useStores();

  const processData = (groupCountArr) => {
    return groupCountArr;
  };

  const handleHoverGroup = (...args) => {
    // Don't fire the action if there's no change
    let hoverGroup = args[1] === null ? null : args[1]['group'];
    if (hoverGroup === configStore.hoverGroup) {
      return;
    }
    configStore.updateHoverGroup(hoverGroup);
  };

  const handleSelected = (...args) => {
    // console.log(args[1], toJS(configStore.selectedGroups));
    const curSelectedGroups = args[1].map((item) => {
      return { group: item.group };
    });
    // Don't fire if the selection is the same
    if (_.isEqual(curSelectedGroups, configStore.selectedGroups)) {
      return;
    } else {
      configStore.updateSelectedGroups(curSelectedGroups);
    }
  };

  const [state, setState] = useState({
    xRange: [0, 100],
    data: {
      table: processData(dataStore.groupCountArr),
      selected: JSON.parse(JSON.stringify(configStore.selectedGroups)),
    },
    signalListeners: {
      hoverGroup: _.throttle(handleHoverGroup, 100),
    },
    dataListeners: {
      selected: handleSelected,
    },
  });

  useEffect(() => {
    // Apply xRange
    let xRange;
    if (configStore.dnaOrAa === DNA_OR_AA.DNA) {
      const coordRanges = toJS(configStore.getCoordinateRanges());
      xRange = [
        coordRanges.reduce((memo, rng) => Math.min(...rng, memo), 30000),
        coordRanges.reduce((memo, rng) => Math.max(...rng, memo), 0),
      ];
    } else if (configStore.dnaOrAa === DNA_OR_AA.AA) {
      // Get the extent of the selected gene/protein
      let residueRanges;
      if (configStore.coordinateMode === COORDINATE_MODES.COORD_GENE) {
        residueRanges = configStore.selectedGene.ranges;
      } else if (
        configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN
      ) {
        residueRanges = configStore.selectedProtein.ranges;
      }
      // Convert NT indices to AA residue indices
      const startNTInd = residueRanges.reduce(
        (memo, rng) => Math.min(...rng, memo),
        30000
      );
      residueRanges = residueRanges.map((rng) =>
        rng.map((ind) => (ind - startNTInd) / 3)
      );
      xRange = [
        residueRanges.reduce((memo, rng) => Math.min(...rng, memo), 30000),
        residueRanges.reduce((memo, rng) => Math.max(...rng, memo), 0),
      ];
    }

    setState({
      ...state,
      xRange,
      data: {
        ...state.data,
        table: processData(dataStore.groupCountArr),
      },
    });
  }, [dataStore.selectedRowsAndDateHash]);

  // Update internal selected groups copy
  useEffect(() => {
    setState({
      ...state,
      data: {
        ...state.data,
        selected: JSON.parse(JSON.stringify(configStore.selectedGroups)),
      },
    });
  }, [configStore.selectedGroups]);

  // Generate x-axis title
  let xLabel = 'WIV04';
  if (configStore.coordinateMode === COORDINATE_MODES.COORD_GENE) {
    xLabel += ', ' + configStore.selectedGene.gene + ' Gene';
  } else if (configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
    xLabel += ', ' + configStore.selectedProtein.protein + ' Protein';
  }

  if (configStore.dnaOrAa === DNA_OR_AA.DNA) {
    xLabel += ' (NT)';
  } else if (configStore.dnaOrAa === DNA_OR_AA.AA) {
    xLabel += ' (AA residue)';
  }

  return (
    <PlotContainer>
      <VegaEmbed
        spec={initialSpec}
        data={state.data}
        width={width}
        signals={{
          xLabel,
          xRange: state.xRange,
          hoverGroup: { group: configStore.hoverGroup },
          posField: configStore.dnaOrAa === DNA_OR_AA.DNA ? 0 : 1,
        }}
        signalListeners={state.signalListeners}
        dataListeners={state.dataListeners}
      />
    </PlotContainer>
  );
});
EntropyPlot.propTypes = {
  width: PropTypes.number,
};
EntropyPlot.defaultProps = {
  width: 100,
};

export default EntropyPlot;
