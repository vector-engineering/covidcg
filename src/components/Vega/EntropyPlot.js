import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { toJS } from 'mobx';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

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

import { geneMap, proteinMap } from '../../utils/gene_protein';
import { throttle } from '../../utils/func';

const PlotContainer = styled.div``;

const EntropyPlot = observer(({ width }) => {
  const vegaRef = useRef();
  const { configStore, dataStore, UIStore, snpDataStore } = useStores();

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

  const filterMap = (map, keyToRemove) => {
    // Used to turn proteinMap and geneMap into an array of objects
    // that can be used with the domain plot
    return Object.keys(map)
      .filter((key) => key !== keyToRemove)
      .reduce((arr, key) => {
        map[key]['ranges'] = map[key]['segments'];
        if (map[key]['ranges'].length > 1) {
          // Collapse ORF1ab ranges
          map[key]['ranges'][0][1] = map[key]['ranges'][1][1];
          map[key]['ranges'].slice(1, 1);
        }
        arr.push(map[key]);
        return arr;
      }, []);
  };

  const processData = () => {
    let snvCounts = toJS(dataStore.groupCounts);
    // Input data from dataStore.groupCounts is in the form
    // [{ group_id: snv_id, counts: int }]
    // Before we pass this into Vega, we also need:
    // 1) The color of each SNV
    // 2) The human-readable name of each SNV
    // 3) The position of the SNV
    // console.log('ENTROPY PROCESS DATA');

    return snvCounts
      .map((record) => {
        let snv = snpDataStore.intToSnv(
          configStore.dnaOrAa,
          configStore.coordinateMode,
          record.group_id
        );

        record.snv = snv.snp_str;
        record.color = snpDataStore.getSnvColor(snv.snp_str);
        record.snvName = snv.name;
        record.pos = snv.pos;
        return record;
      })
      .filter((record) => {
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

  const getDomainPlotHeight = () => {
    // Domain Plot height is calculated as the number of rows times a constant
    let heightConst = 50;
    // There will always be at least 1 row (nullDomain displays when no rows)
    let numRows = 1;

    let geneProteinObj = null;

    // Logic for Primer track
    if (configStore.coordinateMode === COORDINATE_MODES.COORD_PRIMER) {
      if (configStore.selectedPrimers.length > 0) {
        const primerObj = configStore.selectedPrimers;
        let oldInst = primerObj[0].Institution;
        primerObj.forEach((primer, i) => {
          if (primer.Institution !== oldInst) {
            numRows += 1;
            oldInst = primer.Institution;
          }
        });
      }

      if (numRows > 1) heightConst = 15;
      return numRows * heightConst;
    }

    // Logic for Gene/Protein track
    if (configStore.residueCoordinates.length === 0) {
      if (configStore.coordinateMode === COORDINATE_MODES.COORD_GENE) {
        geneProteinObj = filterMap(geneMap, 'All Genes');
      } else if (
        configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN
      ) {
        geneProteinObj = filterMap(proteinMap, 'All Proteins');
      }

      // Greedily get the max number of rows
      geneProteinObj.forEach((geneProtein) => {
        // geneProtein[row] is zero-indexed so add 1 to get total number of rows
        if (geneProtein['row'] + 1 > numRows) {
          numRows = geneProtein['row'] + 1;
        }
      });

      return numRows * heightConst;
    } else if (configStore.coordinateMode === COORDINATE_MODES.COORD_GENE) {
      geneProteinObj = configStore.selectedGene;
    } else if (configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
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

    if (numRows > 1) heightConst = 15;

    return numRows * heightConst;
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

  const getDomains = () => {
    // Apply domains
    const xRange = getXRange();
    const nullDomain = [
      {
        name: 'No Domains Available',
        ranges: [xRange],
        row: 0,
      },
    ];

    if (configStore.coordinateMode === COORDINATE_MODES.COORD_PRIMER) {
      if (configStore.selectedPrimers.length > 0) {
        const selectedPrimers = configStore.selectedPrimers;
        selectedPrimers[0].row = 0;
        let oldInst = selectedPrimers[0].Institution;
        let curRow = 0;
        selectedPrimers.forEach((primer) => {
          if (primer.Institution !== oldInst) {
            curRow += 1;
            oldInst = primer.Institution;
          }
          primer.row = curRow;
          primer.ranges = [[primer.Start, primer.End]];
          primer.name = primer.Name;
        });
        return selectedPrimers;
      } else {
        nullDomain.name = 'No Primers Selected';
        return nullDomain;
      }
    }

    if (configStore.residueCoordinates.length === 0) {
      if (configStore.coordinateMode === COORDINATE_MODES.COORD_GENE) {
        return filterMap(geneMap, 'All Genes');
      } else if (
        configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN
      ) {
        return filterMap(proteinMap, 'All Proteins');
      }
    } else if (configStore.coordinateMode === COORDINATE_MODES.COORD_GENE) {
      return configStore.selectedGene.domains.length > 0
        ? configStore.selectedGene.domains
        : nullDomain;
    } else if (configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
      return configStore.selectedProtein.domains.length > 0
        ? configStore.selectedProtein.domains
        : nullDomain;
    }
  };

  const [state, setState] = useState({
    showWarning: true,
    xRange: getXRange(),
    hoverGroup: null,
    data: { domains: getDomains() },
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

    // console.log(getDomains());
    // console.log(getXRange());
    // console.log(processData());

    setState({
      ...state,
      xRange: getXRange(),
      domainPlotHeight: getDomainPlotHeight(),
      data: {
        ...state.data,
        domains: getDomains(),
        table: processData(),
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
  if (dataStore.numSequencesAfterAllFiltering === 0) {
    return (
      <EmptyPlot height={150}>
        <p>No sequences selected</p>
      </EmptyPlot>
    );
  }

  const entropyPlotHeight = 120;
  const padding = 40;

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
        height={entropyPlotHeight + padding + state.domainPlotHeight}
        signals={{
          totalSequences: dataStore.numSequencesAfterAllFiltering,
          xLabel,
          xRange: state.xRange,
          hoverGroup: state.hoverGroup,
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
