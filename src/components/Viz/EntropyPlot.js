import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { toJS } from 'mobx';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { config } from '../../config';
import { throttle } from '../../utils/func';
import { getMaxReferenceLength } from '../../utils/reference';

import VegaEmbed from '../../react_vega/VegaEmbed';
import WarningBox from '../Common/WarningBox';
import EmptyPlot from '../Common/EmptyPlot';
import SkeletonElement from '../Common/SkeletonElement';
import DropdownButton from '../Buttons/DropdownButton';
import QuestionButton from '../Buttons/QuestionButton';
import {
  PlotOptions,
  OptionSelectContainer,
  OptionInputContainer,
} from './Plot.styles';

import initialSpec from '../../vega_specs/entropy.vg.json';
import {
  ASYNC_STATES,
  COORDINATE_MODES,
  DNA_OR_AA,
  PLOT_DOWNLOAD_OPTIONS,
  GROUPS,
  NORM_MODES,
} from '../../constants/defs.json';
import ExternalLink from '../Common/ExternalLink';

const PlotContainer = styled.div``;

const EntropyPlot = observer(({ width }) => {
  const vegaRef = useRef();
  const {
    configStore,
    dataStore,
    UIStore,
    mutationDataStore,
    plotSettingsStore,
  } = useStores();

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
      dataStore.downloadMutationFrequencies();
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

  const processData = () => {
    let mutationCounts = toJS(dataStore.groupCounts);

    return mutationCounts.filter((record) => {
      return (
        record.group !== GROUPS.OTHER_GROUP &&
        record.group !== GROUPS.REFERENCE_GROUP
      );
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
    configStore.updateSelectedGroups(curSelectedGroups);
  };

  const onChangeEntropyYMode = (event) => {
    const entropyYMode = event.target.value;
    // Default powers
    let entropyYPow = 1.0;
    if (entropyYMode === NORM_MODES.NORM_COUNTS) {
      entropyYPow = 0.5;
    }

    plotSettingsStore.applyPendingChanges({ entropyYMode, entropyYPow });
  };
  const onChangeEntropyYPow = (event) => {
    plotSettingsStore.applyPendingChanges({ entropyYPow: event.target.value });
  };

  // Domain Plot height is calculated as the number of rows times a constant
  const domainPlotRowHeight = 15;
  const getDomainPlotHeight = () => {
    // There will always be at least 1 row (nullDomain displays when no rows)
    let numRows = 1;

    // Logic for Primer track
    if (configStore.coordinateMode === COORDINATE_MODES.COORD_PRIMER) {
      if (configStore.selectedPrimers.length > 0) {
        const primerObj = configStore.selectedPrimers;
        primerObj.forEach((primer) => {
          if (primer.row + 1 > numRows) {
            numRows = primer.row + 1;
          }
        });

        return numRows * domainPlotRowHeight;
      } else {
        return domainPlotRowHeight;
      }
    } else if (
      configStore.coordinateMode === COORDINATE_MODES.COORD_GENE ||
      configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN
    ) {
      // Logic for Gene/Protein track
      let geneProteinObj = null;
      if (configStore.coordinateMode === COORDINATE_MODES.COORD_GENE) {
        geneProteinObj = configStore.selectedGene;
      } else if (
        configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN
      ) {
        geneProteinObj = configStore.selectedProtein;
      }

      // Greedily get the number of rows
      if (geneProteinObj && geneProteinObj.domains.length > 0) {
        geneProteinObj.domains.forEach((domain) => {
          // geneProtein[row] is zero-indexed so add 1 to get total number of rows
          if (domain['row'] + 1 > numRows) {
            numRows = domain['row'] + 1;
          }
        });
      }
    }

    return numRows * domainPlotRowHeight;
  };

  const getXRange = () => {
    // Apply xRange
    let xRange;
    if (configStore.residueCoordinates.length === 0) {
      // If the residue coordinates are empty, then either "All Genes" or
      // "All Proteins" is selected -- so show everything
      xRange = [1, getMaxReferenceLength()];
    } else if (configStore.dnaOrAa === DNA_OR_AA.DNA) {
      const coordRanges = toJS(configStore.getCoordinateRanges());
      // Get the extent of the selected gene/protein
      xRange = [
        coordRanges.reduce(
          (memo, rng) => Math.min(rng[1], rng[2], memo),
          getMaxReferenceLength()
        ),
        coordRanges.reduce((memo, rng) => Math.max(rng[1], rng[2], memo), 0),
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

      xRange = [minResidueIndex, maxResidueIndex];
    }
    return xRange;
  };

  const getDomains = () => {
    // Apply domains
    const xRange = getXRange();
    const nullDomain = [
      {
        name: 'No Domains Available',
        abbr: 'No Domains Available',
        ranges: [xRange],
        nt_ranges: [xRange],
        row: 0,
      },
    ];

    if (configStore.coordinateMode === COORDINATE_MODES.COORD_GENE) {
      return configStore.selectedGene.domains.length > 0
        ? configStore.selectedGene.domains
        : nullDomain;
    } else if (configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
      return configStore.selectedProtein.domains.length > 0
        ? configStore.selectedProtein.domains
        : nullDomain;
    } else {
      return nullDomain;
    }
  };

  const checkIfPrimersOverlap = (primer1, primer2) => {
    return (
      (primer1.Start < primer2.End && primer1.End > primer2.End) ||
      (primer1.Start < primer2.Start && primer1.End > primer2.Start)
    );
  };

  const getPrimers = () => {
    const selectedPrimers = configStore.selectedPrimers;

    if (selectedPrimers.length) {
      selectedPrimers[0].row = 0;
      for (let i = 0; i < selectedPrimers.length; i++) {
        let overlaps = true;
        let curRow = 0;
        const primerToPlace = selectedPrimers[i];

        while (overlaps) {
          const primersInRow = selectedPrimers.filter(
            (primer) => primer.row && primer.row === curRow
          );

          if (primersInRow.length) {
            for (const primer of primersInRow) {
              overlaps = checkIfPrimersOverlap(primer, primerToPlace);
              if (!overlaps) break;
            }
          } else {
            overlaps = false;
          }

          if (overlaps) curRow += 1;
        }

        primerToPlace.row = curRow;
        primerToPlace.ranges = [[primerToPlace.Start, primerToPlace.End]];
        primerToPlace.name = primerToPlace.Name;
        selectedPrimers[i] = primerToPlace;
      }
      return selectedPrimers;
    } else {
      const nullDomain = [
        {
          Institution: 'None',
          Name: 'No Primers Selected',
          ranges: [[0, getMaxReferenceLength()]],
          row: 0,
          Start: 0,
          End: getMaxReferenceLength(),
        },
      ];
      configStore.selectedPrimers = nullDomain;
      return nullDomain;
    }
  };

  const domainsToShow = () => {
    return configStore.coordinateMode === COORDINATE_MODES.COORD_PRIMER
      ? getPrimers()
      : getDomains();
  };

  const [state, setState] = useState({
    showWarning: true,
    xRange: getXRange(),
    hoverGroup: null,
    data: { domains: domainsToShow() },
    domainPlotHeight: getDomainPlotHeight(),
    signalListeners: {
      hoverGroup: throttle(handleHoverGroup, 100),
    },
    dataListeners: {
      selected: handleSelected,
    },
  });

  useEffect(() => {
    setState({
      ...state,
      hoverGroup: { group: configStore.hoverGroup },
    });
  }, [configStore.hoverGroup]);

  // Update internal selected groups copy
  useEffect(() => {
    setState({
      ...state,
      data: {
        ...state.data,
        selected: toJS(configStore.selectedGroups),
      },
    });
  }, [configStore.selectedGroups]);

  const refreshData = () => {
    if (UIStore.caseDataState !== ASYNC_STATES.SUCCEEDED) {
      return;
    }

    // Slightly modify coverage data before we feed it into Vega
    // Since the last data point in the coverage array is set to 0,
    // this will create a glitch in the coverage plot that will show coverage
    // dropping off right at the end
    // Fix this by setting the last coverage data point to the immediate previous
    // coverage data point
    const coverage = mutationDataStore.coverage.slice();
    if (coverage.length >= 2 && coverage[coverage.length - 1].count === 0) {
      coverage[coverage.length - 1].count = coverage[coverage.length - 2].count;
    }

    setState({
      ...state,
      xRange: getXRange(),
      domainPlotHeight: getDomainPlotHeight(),
      data: {
        ...state.data,
        domains: domainsToShow(),
        table: processData(),
        coverage,
      },
    });
  };

  // Refresh data on mount (i.e., tab change) or when data state changes
  useEffect(refreshData, [UIStore.caseDataState]);
  useEffect(refreshData, []);

  // Generate x-axis title
  let xLabel = '';
  if (configStore.dnaOrAa === DNA_OR_AA.DNA) {
    xLabel += 'NT';
  } else if (configStore.dnaOrAa === DNA_OR_AA.AA) {
    xLabel += 'AA residue';
  }
  xLabel += ' (' + configStore.selectedReference;
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
  if (dataStore.numSequencesAfterAllFiltering === 0) {
    return (
      <EmptyPlot height={150}>
        <p>No sequences selected</p>
      </EmptyPlot>
    );
  }

  const entropyPlotHeight = 120;
  const coveragePlotHeight = 40;
  const padding = 40;

  return (
    <PlotContainer>
      <WarningBox show={state.showWarning} onDismiss={onDismissWarning}>
        {config.site_title} plots reflect data contributed to{' '}
        {config.data_provider} and are therefore impacted by the sequence
        coverage in each country.
        {config['virus'] === 'sars2' && (
          <>
            {' '}
            For example, systematic errors are sometimes observed specific to
            particular labs or methods (
            <ExternalLink href="https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473/14">
              https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473/14
            </ExternalLink>
            ,{' '}
            <ExternalLink href="https://doi.org/10.1371/journal.pgen.1009175">
              https://doi.org/10.1371/journal.pgen.1009175
            </ExternalLink>
            ). Users are advised to consider these errors in their high
            resolution analyses.
          </>
        )}
      </WarningBox>
      <PlotOptions>
        <OptionSelectContainer>
          <label>
            Show:
            <select
              value={plotSettingsStore.entropyYMode}
              onChange={onChangeEntropyYMode}
            >
              <option value={NORM_MODES.NORM_COUNTS}>Counts</option>
              <option value={NORM_MODES.NORM_PERCENTAGES}>Percents</option>
              <option value={NORM_MODES.NORM_COVERAGE_ADJUSTED}>
                Percents (coverage adjusted)
              </option>
            </select>
          </label>
          <QuestionButton
            data-tip={`
          <ul>
            <li>
              Counts: show raw mutation counts
            </li>
            <li>
              Percents: show mutation counts as a percentage of 
              the total number of sequences selected
            </li>
            <li>
              Percents (coverage adjusted): show mutation counts as a 
              percentage of sequences with coverage at the mutation position
            </li>
          </ul>
          `}
            data-html="true"
            data-for="main-tooltip"
          />
        </OptionSelectContainer>

        <OptionInputContainer style={{ marginLeft: '10px' }}>
          <label>
            Y-scale:
            <input
              value={plotSettingsStore.entropyYPow}
              type="number"
              step={0.1}
              min={0}
              onChange={onChangeEntropyYPow}
            ></input>
          </label>
          <QuestionButton
            data-tip={`
          <ul>
            <li>
              Y-scale: adjust the power of the y-axis. 
            </li>
            <li>
              For Y-scale = 1, the y-axis scale is linear.
            </li>
            <li>
              For Y-scale < 1, lower values are more visible.
            </li>
            <li>
              For Y-scale > 1, lower values are less visible.
            </li>
          </ul>
          `}
            data-html="true"
            data-for="main-tooltip"
          />
        </OptionInputContainer>
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
        height={
          entropyPlotHeight +
          padding +
          state.domainPlotHeight +
          padding +
          coveragePlotHeight
        }
        signals={{
          yMode: plotSettingsStore.entropyYMode,
          yScaleExponent: plotSettingsStore.entropyYPow,
          totalSequences: dataStore.numSequencesAfterAllFiltering,
          xLabel,
          xRange: state.xRange,
          hoverGroup: state.hoverGroup,
          numDomainRows: state.domainPlotHeight / domainPlotRowHeight,
          domainPlotHeight: state.domainPlotHeight,
          posField:
            configStore.dnaOrAa === DNA_OR_AA.DNA &&
            configStore.residueCoordinates.length !== 0 &&
            configStore.coordinateMode !== COORDINATE_MODES.COORD_PRIMER
              ? 0
              : 1,
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
  height: PropTypes.number,
};
EntropyPlot.defaultProps = {
  width: 100,
  height: 220,
};

export default EntropyPlot;
