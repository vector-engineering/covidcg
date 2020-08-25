import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { toJS } from 'mobx';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import _ from 'underscore';

import VegaEmbed from '../../react_vega/VegaEmbed';
import EmptyPlot from '../Common/EmptyPlot';
import SkeletonElement from '../Common/SkeletonElement';
import DropdownButton from '../Buttons/DropdownButton';
import { PlotTitle, PlotOptions } from './Plot.styles';

import initialSpec from '../../vega_specs/entropy.vg.json';
import { ASYNC_STATES } from '../../constants/UI';
import { COORDINATE_MODES, DNA_OR_AA } from '../../constants/config';
import { PLOT_DOWNLOAD_OPTIONS } from '../../constants/download';

const PlotContainer = styled.div``;

const EntropyPlot = observer(({ width }) => {
  const vegaRef = useRef();
  const { configStore, dataStore, UIStore } = useStores();

  const handleDownloadSelect = (option) => {
    // console.log(option);
    // TODO: use the plot options and configStore options to build a more descriptive filename
    //       something like new_lineages_by_day_S_2020-05-03-2020-05-15_NYC.png...
    if (option === PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_DATA) {
      dataStore.downloadSnvFrequencies();
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

  const processData = (snvCounts) => {
    return snvCounts;
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

  const getXRange = () => {
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
    return xRange;
  };

  const [state, setState] = useState({
    xRange: getXRange(),
    data: {
      table: processData(toJS(dataStore.groupCountDateFilteredArr)),
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
    if (UIStore.aggCaseDataState !== ASYNC_STATES.SUCCEEDED) {
      return;
    }

    setState({
      ...state,
      xRange: getXRange(),
      data: {
        ...state.data,
        table: processData(toJS(dataStore.groupCountDateFilteredArr)),
      },
    });
  }, [UIStore.aggCaseDataState]);

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
  let xLabel = '';
  if (configStore.dnaOrAa === DNA_OR_AA.DNA) {
    xLabel += 'NT';
  } else if (configStore.dnaOrAa === DNA_OR_AA.AA) {
    xLabel += 'AA residue';
  }
  xLabel += ' (WIV04';
  if (configStore.coordinateMode === COORDINATE_MODES.COORD_GENE) {
    xLabel += ', ' + configStore.selectedGene.gene + ' Gene';
  } else if (configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
    xLabel += ', ' + configStore.selectedProtein.protein + ' Protein';
  }
  xLabel += ')';

  if (
    UIStore.caseDataState === ASYNC_STATES.STARTED ||
    UIStore.aggCaseDataState === ASYNC_STATES.STARTED
  ) {
    return (
      <div
        style={{
          paddingTop: '12px',
          paddingRight: '24px',
          paddingLeft: '12px',
          paddingBottom: '24px',
        }}
      >
        <SkeletonElement delay={2} height={150} />
      </div>
    );
  }

  // If we have no rows, then return an empty element
  // We'll always have the "reference" row, so no rows = 1 row
  if (dataStore.selectedAccessionIds.length === 0) {
    return (
      <EmptyPlot height={150}>
        <p>No sequences selected</p>
      </EmptyPlot>
    );
  }

  return (
    <PlotContainer>
      <PlotOptions>
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
        data={state.data}
        width={width}
        signals={{
          totalSequences: dataStore.selectedAccessionIds.length,
          xLabel,
          xRange: state.xRange,
          hoverGroup: { group: configStore.hoverGroup },
          posField: configStore.dnaOrAa === DNA_OR_AA.DNA ? 0 : 1,
        }}
        signalListeners={state.signalListeners}
        dataListeners={state.dataListeners}
        actions={false}
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
