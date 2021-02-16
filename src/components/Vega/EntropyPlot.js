import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { toJS } from 'mobx';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import _ from 'underscore';

import VegaEmbed from '../../react_vega/VegaEmbed';
import WarningBox from '../Common/WarningBox';
import EmptyPlot from '../Common/EmptyPlot';
import SkeletonElement from '../Common/SkeletonElement';
import DropdownButton from '../Buttons/DropdownButton';
import { PlotOptions } from './Plot.styles';

import initialSpec from '../../vega_specs/entropy.vg.json';
import {
  ASYNC_STATES,
  COORDINATE_MODES,
  DNA_OR_AA,
  PLOT_DOWNLOAD_OPTIONS,
  GROUPS,
} from '../../constants/defs.json';
import ExternalLink from '../Common/ExternalLink';

const PlotContainer = styled.div``;

const EntropyPlot = observer(({ width }) => {
  const vegaRef = useRef();
  const { configStore, dataStore, UIStore, plotSettingsStore } = useStores();

  const onDismissWarning = () => {
    setState({
      ...state,
      showWarning: false,
    });
  };

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
    return snvCounts.filter((group) => {
      return group[0] !== GROUPS.OTHER_GROUP;
    });
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
    if (configStore.residueCoordinates.length === 0) {
      // If the residue coordinates are empty, then either "All Genes" or
      // "All Proteins" is selected -- so show everything
      xRange = [1, 30000];
    } else if (configStore.dnaOrAa === DNA_OR_AA.DNA) {
      const coordRanges = toJS(configStore.getCoordinateRanges());
      xRange = [
        coordRanges.reduce((memo, rng) => Math.min(...rng, memo), 30000),
        coordRanges.reduce((memo, rng) => Math.max(...rng, memo), 0),
      ];
    } else if (configStore.dnaOrAa === DNA_OR_AA.AA) {
      // Get the extent of the selected gene/protein
      let geneOrProteinObj;
      let residueCoordinates = toJS(configStore.residueCoordinates);
      if (configStore.coordinateMode === COORDINATE_MODES.COORD_GENE) {
        geneOrProteinObj = configStore.selectedGene;
      } else if (
        configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN
      ) {
        geneOrProteinObj = configStore.selectedProtein;
      }

      // Find the smallest and largest residue index
      // from the selected residue coordinates
      const minResidueIndex = residueCoordinates.reduce(
        (minIndex, rng) => Math.min(...rng, minIndex),
        geneOrProteinObj.len_aa
      );
      const maxResidueIndex = residueCoordinates.reduce(
        (minIndex, rng) => Math.max(...rng, minIndex),
        1
      );

      if (configStore.dnaOrAa === DNA_OR_AA.DNA) {
        // Find the AA range that the minimum and maximum AA index fall in,
        // and then get the NT coordinates (from the start of the codon)
        let startNTInd, endNTInd;
        geneOrProteinObj.aa_ranges.forEach((aaRange, ind) => {
          if (minResidueIndex >= aaRange[0] && minResidueIndex <= aaRange[1]) {
            // Get the matching NT range, add residues * 3
            startNTInd =
              geneOrProteinObj.ranges[ind][0] +
              (minResidueIndex - aaRange[0]) * 3;
          }
          if (maxResidueIndex >= aaRange[0] && maxResidueIndex <= aaRange[1]) {
            // Get the matching NT range, add residues * 3 (to end of codon)
            endNTInd =
              geneOrProteinObj.ranges[ind][0] +
              2 +
              (maxResidueIndex - aaRange[0]) * 3;
          }
        });
        xRange = [startNTInd, endNTInd];
      } else if (configStore.dnaOrAa === DNA_OR_AA.AA) {
        xRange = [minResidueIndex, maxResidueIndex];
      }
    }
    return xRange;
  };

  const [state, setState] = useState({
    showWarning: true,
    xRange: getXRange(),
    data: {
      table: processData(toJS(dataStore.countsPerGroupDateFiltered)),
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
    if (UIStore.caseDataState !== ASYNC_STATES.SUCCEEDED) {
      return;
    }

    setState({
      ...state,
      xRange: getXRange(),
      data: {
        ...state.data,
        table: processData(toJS(dataStore.countsPerGroupDateFiltered)),
      },
    });
  }, [UIStore.caseDataState, plotSettingsStore.entropyMinCount]);

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
    xLabel += ', ' + configStore.selectedGene.name + ' Gene';
  } else if (configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
    xLabel += ', ' + configStore.selectedProtein.name + ' Protein';
  }
  xLabel += ')';

  if (UIStore.caseDataState === ASYNC_STATES.STARTED) {
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
  if (dataStore.filteredCaseData.length === 0) {
    return (
      <EmptyPlot height={150}>
        <p>No sequences selected</p>
      </EmptyPlot>
    );
  }

  return (
    <PlotContainer>
      <WarningBox show={state.showWarning} onDismiss={onDismissWarning}>
        Systematic errors are sometimes observed specific to particular labs or
        methods (
        <ExternalLink href="https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473/14">
          https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473/14
        </ExternalLink>
        ,{' '}
        <ExternalLink href="https://doi.org/10.1371/journal.pgen.1009175">
          https://doi.org/10.1371/journal.pgen.1009175
        </ExternalLink>
        ), users are advised to consider these errors in their high resolution
        analyses.
      </WarningBox>
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
          totalSequences: dataStore.filteredCaseData.length,
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
