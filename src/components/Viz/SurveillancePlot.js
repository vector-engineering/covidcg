import React, { useRef, useEffect, useState } from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { format } from 'd3-format';
import { config } from '../../config';

import {
  SORT_DIRECTIONS,
  PLOT_DOWNLOAD_OPTIONS,
  SURV_GRAPH_MODE,
} from '../../constants/defs.json';

import ReactTooltip from 'react-tooltip';
import VegaEmbed from '../../react_vega/VegaEmbed';
import {
  PlotOptions,
  OptionSelectContainer,
  OptionInputContainer,
} from './Plot.styles';
import DropdownButton from '../Buttons/DropdownButton';
import QuestionButton from '../Buttons/QuestionButton';
import WarningBox from '../Common/WarningBox';

import initialSpec from '../../vega_specs/surveillance.vg.json';
import stackedSpec from '../../vega_specs/surveillance_stacked.vg.json';

import {
  PlotContainer,
  PlotTitle,
  HelpText,
  OptionsColumn,
  OptionsRow,
  PlotAndLegend,
  Legend,
  LegendTitle,
  LegendList,
  LegendItemContainer,
  LegendItemName,
  LegendItemCounts,
  CollapseButton,
} from './SurveillancePlot.styles';

const countsFormatter = format('.2s');

const LegendItem = ({ legendItem, legendHover, setLegendHover }) => {
  const onMouseEnter = () => {
    // console.log('enter', legendItem.group);
    setLegendHover(legendItem.group);
  };

  const onMouseOut = () => {
    // console.log('leave', legendItem.group);
    setLegendHover(null);
  };

  let groupText = legendItem.group;

  return (
    <LegendItemContainer
      hovered={legendHover.includes(legendItem.group)}
      barColor={legendItem.color}
      onMouseEnter={onMouseEnter}
      onMouseLeave={onMouseOut}
    >
      <LegendItemName>{groupText}</LegendItemName>
      <LegendItemCounts>{countsFormatter(legendItem.counts)}</LegendItemCounts>
    </LegendItemContainer>
  );
};
LegendItem.propTypes = {
  legendItem: PropTypes.object,
  legendHover: PropTypes.arrayOf(PropTypes.string),
  setLegendHover: PropTypes.func,
};

const SurveillanceLegend = observer(
  ({ legendItems, legendHover, setLegendHover }) => {
    const { plotSettingsStore } = useStores();
    const legendElements = [];
    legendItems.forEach((legendItem) => {
      legendElements.push(
        <LegendItem
          key={`surv-legend-${legendItem.group}`}
          legendItem={legendItem}
          legendHover={legendHover}
          setLegendHover={setLegendHover}
        />
      );
    });

    let legendTitleText = plotSettingsStore.surveillanceMode;
    // Capitalize first letter
    legendTitleText =
      legendTitleText.charAt(0).toUpperCase() + legendTitleText.slice(1);

    return (
      <Legend>
        <LegendTitle>{legendTitleText}</LegendTitle>
        <LegendList>{legendElements}</LegendList>
      </Legend>
    );
  }
);
SurveillanceLegend.propTypes = {
  legendItems: PropTypes.arrayOf(PropTypes.object),
  legendHover: PropTypes.arrayOf(PropTypes.string),
  setLegendHover: PropTypes.func,
};

const SurveillanceSettings = observer(
  ({ groupName, vegaRef, refreshValidGroups }) => {
    const { plotSettingsStore } = useStores();

    // Build tooltips once we're mounted
    useEffect(() => {
      ReactTooltip.rebuild();
    }, []);

    const onChangeSortField = (e) => {
      const oldSortField = plotSettingsStore.surveillanceSortField;
      const newSortField = e.target.value;

      plotSettingsStore.applyPendingChanges({
        surveillanceSortField: newSortField,
      });

      // Manually trigger a refresh on the legend data
      if (oldSortField !== newSortField) {
        vegaRef.current.runWhenComplete(refreshValidGroups);
      }
    };
    const onChangeSortDirection = (e) => {
      const oldSortDirection = plotSettingsStore.surveillanceSortDirection;
      const newSortDirection = e.target.value;

      plotSettingsStore.applyPendingChanges({
        surveillanceSortDirection: e.target.value,
      });

      // Manually trigger a refresh on the legend data
      // console.log(oldSortDirection, newSortDirection);
      if (oldSortDirection !== newSortDirection) {
        vegaRef.current.runWhenComplete(refreshValidGroups);
      }
    };
    const onChangeDisplayMinCounts = (e) => {
      plotSettingsStore.applyPendingChanges({
        surveillanceDisplayMinCounts: e.target.value,
      });
    };
    const onChangeDisplayMinPercent = (e) => {
      plotSettingsStore.applyPendingChanges({
        surveillanceDisplayMinPercent: e.target.value,
      });
    };
    const onChangeSigMinCounts = (e) => {
      plotSettingsStore.applyPendingChanges({
        surveillanceSigMinCounts: e.target.value,
      });
    };
    const onChangeSigMinPercent = (e) => {
      plotSettingsStore.applyPendingChanges({
        surveillanceSigMinPercent: e.target.value,
      });
    };
    const onChangeSigMinR = (e) => {
      plotSettingsStore.applyPendingChanges({
        surveillanceSigMinR: e.target.value,
      });
    };

    return (
      <div>
        <OptionsRow>
          Sort by{' '}
          <OptionSelectContainer>
            <label>
              <select
                value={plotSettingsStore.surveillanceSortField}
                onChange={onChangeSortField}
              >
                <option value={'group'}>Name</option>
                <option value={'counts'}>Counts</option>
              </select>
            </label>
          </OptionSelectContainer>
          <OptionSelectContainer>
            <label>
              <select
                value={plotSettingsStore.surveillanceSortDirection}
                onChange={onChangeSortDirection}
              >
                <option value={SORT_DIRECTIONS.SORT_ASC}>Ascending</option>
                <option value={SORT_DIRECTIONS.SORT_DESC}>Descending</option>
              </select>
            </label>
          </OptionSelectContainer>
        </OptionsRow>
        <OptionsRow>
          <span>Displayed {groupName}</span>
          <QuestionButton
            data-tip={`<p>To improve readability, ${groupName} which do not meet these conditions in any of the six continents are removed from the plots. For example, if one ${groupName} fails these conditions in 5 continents but passes in the last one, it will still be shown. If one ${groupName} fails these conditions in all 6 continents, it will be excluded.</p><p>Min Counts is defined as the minimum number of sequences with this ${groupName} in a given continent on a given week</p><p>Min Percent is defined as the minimum percent of sequences with this ${groupName} in a given continent on a given week</p>`}
            data-html="true"
            data-place="right"
            data-for="main-tooltip"
            style={{ marginRight: 8 }}
          />
          <OptionInputContainer>
            <label>
              Min Counts
              <input
                type="number"
                value={plotSettingsStore.surveillanceDisplayMinCounts}
                onChange={onChangeDisplayMinCounts}
                min={0}
                step={1}
              />
            </label>
          </OptionInputContainer>
          <OptionInputContainer>
            <label>
              Min Percent
              <input
                type="number"
                value={plotSettingsStore.surveillanceDisplayMinPercent}
                onChange={onChangeDisplayMinPercent}
                min={0}
                max={1}
                step={0.01}
              />
            </label>
          </OptionInputContainer>
        </OptionsRow>
        <OptionsRow>
          <span>Highlighted {groupName}</span>
          <QuestionButton
            data-tip={`<p>Choose ${groupName} to highlight and display in the legend. ${groupName} which pass these conditions in at least one of six continents will be highlighted. Non-highlighted ${groupName} are still displayed but not colored.</p><p>Min Counts is defined as the minimum number of sequences with this ${groupName} in a given continent on a given week</p><p>Min Percent is defined as the minimum percent of sequences with this ${groupName} in a given continent on a given week</p><p>Min Correlation is defined as the minimum Pearson correlation of the % sequences of this ${groupName}. For the correlation calculations, the last ${config.surv_end_date_days_ago} days from today are omitted, to reduce noise in this metric.</p>`}
            data-html="true"
            data-place="right"
            data-for="main-tooltip"
            style={{ marginRight: 8 }}
          />
          <OptionInputContainer>
            <label>
              Min Counts
              <input
                type="number"
                value={plotSettingsStore.surveillanceSigMinCounts}
                onChange={onChangeSigMinCounts}
                min={0}
                step={1}
              />
            </label>
          </OptionInputContainer>
          <OptionInputContainer>
            <label>
              Min Percent
              <input
                type="number"
                value={plotSettingsStore.surveillanceSigMinPercent}
                onChange={onChangeSigMinPercent}
                min={0}
                max={1}
                step={0.01}
              />
            </label>
          </OptionInputContainer>
          <OptionInputContainer>
            <label>
              Min Correlation
              <input
                type="number"
                value={plotSettingsStore.surveillanceSigMinR}
                onChange={onChangeSigMinR}
                min={-1}
                max={1}
                step={0.01}
              />
            </label>
          </OptionInputContainer>
        </OptionsRow>
      </div>
    );
  }
);
SurveillanceSettings.propTypes = {
  groupName: PropTypes.string.isRequired,
  vegaRef: PropTypes.object.isRequired,
  refreshValidGroups: PropTypes.func.isRequired,
};

const SurveillancePlot = observer(({ width }) => {
  const { plotSettingsStore, surveillanceDataStore } = useStores();
  const vegaRef = useRef();
  const hiddenLink = useRef();

  // Build tooltips once we're mounted
  useEffect(() => {
    ReactTooltip.rebuild();
  }, []);

  const handleValidGroups = (...args) => {
    // console.log('valid', args);
    setState({
      ...state,
      legendItems: args[1],
    });
  };

  // Long story for this one
  // We need to pull out the "valid_groups" dataset from the Vega plot
  // in order to render our legend
  // To get this, we use the dataListener to get changes in this dataset
  // so we can update the legend accordingly.
  // HOWEVER, on the initial plot rendering, the data callback doesn't fire
  // and I can't figure out a way to configure the callback to fire after initialization
  // So, the hack around this is to call this function once the plot initializes and
  // finishes it's first data flow. We grab the "valid_groups" dataset and pipe it
  // through the callback function as if the callback function fired in the beginning
  const refreshValidGroups = () => {
    vegaRef.current.getData('valid_groups_color', (validGroups) => {
      handleValidGroups('valid_groups_color', validGroups);
    });
  };

  const handleHoverGroups = (...args) => {
    // console.log('hover', args);

    const legendHover = [];
    args[1].forEach((item) => {
      legendHover.push(item.group);
    });

    plotSettingsStore.applyPendingChanges({
      surveillanceLegendHover: legendHover,
    });
  };

  const [state, setState] = useState({
    dataListeners: {
      valid_groups_color: handleValidGroups,
      tooltip_group: handleHoverGroups,
    },
    data: {
      group_reg: surveillanceDataStore.surv_group_regression,
      group_counts: surveillanceDataStore.surv_group_counts,
      hover_legend: [],
    },
    vegaSpec: stackedSpec,
    signalListeners: {},
    legendItems: [],
  });

  const onToggleShowSettings = () => {
    plotSettingsStore.applyPendingChanges({
      surveillanceShowSettings: !plotSettingsStore.surveillanceShowSettings,
    });
  };

  const onDismissWarning = () => {
    plotSettingsStore.applyPendingChanges({ surveillanceShowWarning: false });
  };

  const handleDownloadSelect = (option) => {
    if (option === PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_DATA) {
      surveillanceDataStore.downloadGroupCounts();
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

  // const onChangeMode = (e) => {
  //   plotSettingsStore.applyPendingChanges({ surveillanceMode: e.target.value });
  // };

  const onChangeGraphMode = (e) => {
    plotSettingsStore.applyPendingChanges({
      surveillanceGraphMode: e.target.value,
    });
  };

  const setLegendHover = (group) => {
    // console.log('legend hover', group);

    let legendHover;
    if (group === null) {
      legendHover = [];
    } else if (plotSettingsStore.surveillanceLegendHover.includes(group)) {
      return;
    } else {
      legendHover = [group];
    }

    plotSettingsStore.applyPendingChanges({
      surveillanceLegendHover: legendHover,
    });

    setState({
      ...state,
      data: {
        hover_legend: legendHover.map((group) => {
          return { group };
        }),
      },
    });
  };

  let groupName = plotSettingsStore.surveillanceMode;
  // Capitalize first letter
  groupName = groupName.charAt(0).toUpperCase() + groupName.slice(1);

  const getXLabelFormat = () => {
    if (config.virus === 'sars2') {
      return '%m-%d';
    } else {
      return '%Y-%m';
    }
  };

  // For snapping to nearest points, on the x (time) axis
  const getTimeSensitivity = () => {
    const day = 24 * 60 * 60 * 1000;
    if (config.surv_period === 'W') {
      return day * 3.5;
    } else if (config.surv_period === 'M') {
      return day * 15;
    } else if (config.surv_period === 'Y') {
      return day * 365;
    } else {
      return day * 100;
    }
  };

  useEffect(() => {
    setState({
      ...state,
      data: {
        group_reg: surveillanceDataStore.surv_group_regression,
        group_counts: surveillanceDataStore.surv_group_counts,
      },
    });
  }, [plotSettingsStore.surveillanceMode]);

  let survPeriodText = '';
  if (config.surv_period === 'W') {
    survPeriodText = 'week';
  } else if (config.surv_period === 'M') {
    survPeriodText = 'month';
  } else if (config.surv_period === 'Y') {
    survPeriodText = 'year';
  }

  let vegaSpec;
  if (plotSettingsStore.surveillanceGraphMode === SURV_GRAPH_MODE.LINE) {
    vegaSpec = initialSpec;
  } else if (
    plotSettingsStore.surveillanceGraphMode === SURV_GRAPH_MODE.STACK
  ) {
    vegaSpec = stackedSpec;
  }

  console.log(state.data);

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
      </PlotOptions>
      <HelpText>
        Only data from{' '}
        {Object.prototype.hasOwnProperty.call(config, 'surv_start_date')
          ? config.surv_start_date
          : `the last ${config.surv_start_date_days_ago} days`}{' '}
        is shown, and sequence counts are grouped by {survPeriodText} to reduce
        noise. Please note that the most recent data (most recent month marked
        by darker-colored band) is sparser due to lags in time between sample
        collection and submission.
      </HelpText>
      <HelpText>
        {groupName}s that do not meet the conditions defined by &quot;Displayed{' '}
        {groupName}&quot; in all six continents are filtered out of this plot.
        &quot;Highlighted {groupName}&quot; meeting the user-defined conditions
        in at least one out of the six continents are shown in the legend to the
        left. Hover over {groupName}s in the legend, or near them in the plots,
        to highlight the {groupName} across all plots.
      </HelpText>
      <WarningBox
        show={plotSettingsStore.surveillanceShowWarning}
        onDismiss={onDismissWarning}
      >
        {config.site_title} plots reflect data contributed to{' '}
        {config.data_provider} and are therefore impacted by the sequence
        coverage in each country. Increased prevalence of any {groupName} does
        not, on its own, suggest an increase in transmissibility.
      </WarningBox>
      {/* Only show these surveillance settings for SAR2 */}
      {config.virus === 'sars2' && (
        <PlotOptions style={{ marginLeft: 20 }}>
          <OptionsColumn>
            <OptionsRow>
              <CollapseButton onClick={onToggleShowSettings}>
                {plotSettingsStore.surveillanceShowSettings ? 'Hide' : 'Show'}{' '}
                Plot Settings <span className="caret"></span>
              </CollapseButton>
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
            </OptionsRow>
            {plotSettingsStore.surveillanceShowSettings && (
              <SurveillanceSettings
                groupName={groupName}
                vegaRef={vegaRef}
                refreshValidGroups={refreshValidGroups}
              />
            )}
          </OptionsColumn>
        </PlotOptions>
      )}
      <PlotOptions>
        <OptionsColumn>
          <OptionSelectContainer>
            <label>
              <select
                value={plotSettingsStore.surveillanceGraphMode}
                onChange={onChangeGraphMode}
              >
                <option value={SURV_GRAPH_MODE.LINE}>Line</option>
                <option value={SURV_GRAPH_MODE.STACK}>Stacked Bar</option>
              </select>
            </label>
          </OptionSelectContainer>
        </OptionsColumn>
      </PlotOptions>
      <PlotAndLegend>
        <SurveillanceLegend
          legendItems={state.legendItems}
          legendHover={plotSettingsStore.surveillanceLegendHover}
          setLegendHover={setLegendHover}
        />
        <VegaEmbed
          ref={vegaRef}
          spec={vegaSpec}
          data={state.data}
          signals={{
            sortField: plotSettingsStore.surveillanceSortField,
            sortDirection:
              plotSettingsStore.surveillanceSortDirection ===
              SORT_DIRECTIONS.SORT_ASC
                ? 'ascending'
                : 'descending',
            display_min_counts: plotSettingsStore.surveillanceDisplayMinCounts,
            display_min_percent:
              plotSettingsStore.surveillanceDisplayMinPercent,
            sig_min_counts: plotSettingsStore.surveillanceSigMinCounts,
            sig_min_percent: plotSettingsStore.surveillanceSigMinPercent,
            sig_min_r: plotSettingsStore.surveillanceSigMinR,
            xLabelFormat: getXLabelFormat(),
            time_sensitivity: getTimeSensitivity(),
          }}
          dataListeners={state.dataListeners}
          width={width - 80}
          actions={false}
          onComplete={refreshValidGroups}
        />
      </PlotAndLegend>
    </PlotContainer>
  );
});
SurveillancePlot.propTypes = {
  width: PropTypes.number,
};
SurveillancePlot.defaultProps = {
  width: 1000,
};

export default SurveillancePlot;
