import React from 'react';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';
import useDimensions from 'react-use-dimensions';
import {
  GROUP_MUTATION,
  DNA_OR_AA,
  COORDINATE_MODES,
  TABS,
} from '../../constants/defs.json';

import KBD from '../Common/KBD';
import TabIndicator from '../Common/TabIndicator';
import SelectBoxText from '../Common/SelectBoxText';
import AccordionWrapper from '../Common/AccordionWrapper';

import VegaStackedBars from '../Viz/GroupStackPlot';
import LocationGroupPlot from '../Viz/LocationGroupPlot';
import EntropyPlot from '../Viz/EntropyPlot';
import CooccurrencePlot from '../Viz/CooccurrencePlot';
import NumSeqPerLocationBar from '../Viz/NumSeqPerLocationBar';
import NumSeqPerLocationLine from '../Viz/NumSeqPerLocationLine';
import MutationStructureViewer from '../Viz/MutationStructureViewer';

const CompareGroupsTabContainer = styled.div`
  padding-top: 10px;
  padding-bottom: 50px;
`;

const CompareGroupsTab = observer(() => {
  const { configStore } = useStores();
  const [ref, { width }] = useDimensions();

  const renderEntropyPlot = () => {
    if (configStore.groupKey !== GROUP_MUTATION) {
      return null;
    }

    return (
      <AccordionWrapper
        title={`${configStore.getGroupLabel()} Frequencies`}
        defaultCollapsed={false}
        maxHeight={'500px'}
        helpText={
          <ul>
            <li>
              This plot shows the frequency of {configStore.getGroupLabel()}s
              along the {configStore.selectedReference} genome (
              {configStore.dnaOrAa === DNA_OR_AA.DNA ? 'NT' : 'AA'} Mode:{' '}
              {configStore.dnaOrAa === DNA_OR_AA.DNA
                ? 'Genomic Coordinates'
                : 'Residue Indices'}
              ).
            </li>
            <li>
              <i>Click</i> to select one, or hold <KBD>Shift</KBD> and{' '}
              <i>click</i> to select multiple {configStore.getGroupLabel()}s.
            </li>
            <li>
              Selected {configStore.getGroupLabel()}s will be highlighted in the
              plots and table below, as well as in the{' '}
              <TabIndicator tab={TABS.TAB_COMPARE_LOCATIONS}>
                Compare Locations
              </TabIndicator>{' '}
              tab.
            </li>
          </ul>
        }
      >
        <EntropyPlot width={width - 200} />
      </AccordionWrapper>
    );
  };

  const renderCooccurrencePlot = () => {
    if (configStore.groupKey !== GROUP_MUTATION) {
      return null;
    }

    return (
      <AccordionWrapper
        title={`${configStore.getGroupLabel()} Co-occurrence`}
        defaultCollapsed={false}
        maxHeight={'600px'}
        helpText={
          <ul>
            <li>
              This plot shows {configStore.getGroupLabel()}s that are
              co-occurring with the selected {configStore.getGroupLabel()}s
              (labels on the y-axis).
            </li>
            <li>
              Mutations are counted per-co-occurrence, i.e., one matching
              sequence may count towards multiple mutations. Only the relative
              counts between co-occurring mutations within the same bar should
              be interpreted – the sum of all mutation counts per bar is not
              meaningful.
            </li>
            <li>
              <i>Click</i> on a mutation bar, or on a mutation y-axis label, to
              select or deselect a mutation. If the mutation is already
              selected, then clicking will deselect it. If it is not already
              selected, then clicking on the mutation will add it to the
              existing selection.
            </li>
            <li>
              Mutation frequencies can be shown as raw counts{' '}
              <SelectBoxText>Counts</SelectBoxText>, or counts normalized
              between each selected {configStore.getGroupLabel()}s or
              combination of {configStore.getGroupLabel()}s{' '}
              <SelectBoxText>Normalized Counts</SelectBoxText>.
            </li>
          </ul>
        }
      >
        <CooccurrencePlot width={width - 300} />
      </AccordionWrapper>
    );
  };

  const renderGroupStackPlot = () => {
    return (
      <AccordionWrapper
        title={`${configStore.getGroupLabel()} Plot`}
        defaultCollapsed={false}
        maxHeight={'1200px'}
        helpText={
          <ul>
            <li>
              The plot shows sequences grouped by their respective{' '}
              <b>{configStore.getGroupLabel()}</b> and plotted over time.
            </li>
            {configStore.groupKey === GROUP_MUTATION && (
              <li>
                In <b>Mutation Mode</b>, sequences are split into two groups:
                those that have the mutation (or combination of mutations), and
                those that don&apos;t.
              </li>
            )}
            <li>
              <i>Click</i> to select one, or hold <KBD>Shift</KBD> and{' '}
              <i>click</i> to select multiple {configStore.getGroupLabel()}s.
            </li>
            <li>
              Selected {configStore.getGroupLabel()}s will be highlighted in the
              legend and table below, as well as in the{' '}
              <TabIndicator tab={TABS.TAB_COMPARE_LOCATIONS}>
                Compare Locations
              </TabIndicator>{' '}
              tab.
            </li>
            <li>
              <i>Click</i> and <i>drag</i> on the lower plot (&quot;All
              Seqs&quot;) to zoom in on a specific date range.
            </li>
            <li>
              {configStore.getGroupLabel()} frequencies can be shown as new
              occurrences <SelectBoxText>New</SelectBoxText>, or as a cumulative
              count over time <SelectBoxText>Cumulative</SelectBoxText>
            </li>
            <li>
              {configStore.getGroupLabel()}s can be plotted as raw counts{' '}
              <SelectBoxText>Counts</SelectBoxText>, or as percentages
              normalized to the total number of sequences at that time
              <SelectBoxText>Percentages</SelectBoxText>
            </li>
            <li>
              Sequences are plotted over time by the sample collection date.
              Sequences can be grouped into time bins by{' '}
              <SelectBoxText>Day</SelectBoxText>,{' '}
              <SelectBoxText>Week</SelectBoxText>, or{' '}
              <SelectBoxText>Month</SelectBoxText>
            </li>
          </ul>
        }
      >
        <VegaStackedBars width={width - 150} />
        {configStore.groupKey !== GROUP_MUTATION && (
          <LocationGroupPlot width={width - 250} />
        )}
      </AccordionWrapper>
    );
  };

  const renderMutationStructureViewer = () => {
    // Only valid in AA mode, gene/protein mode
    if (
      configStore.dnaOrAa !== DNA_OR_AA.AA ||
      (configStore.coordinateMode !== COORDINATE_MODES.COORD_GENE &&
        configStore.coordinateMode !== COORDINATE_MODES.COORD_PROTEIN)
    ) {
      return null;
    }

    return (
      <AccordionWrapper
        title="Mutations projected onto protein structure"
        defaultCollapsed={false}
        maxHeight={'500px'}
        helpText={
          <ul>
            <li>...</li>
          </ul>
        }
      >
        <MutationStructureViewer />
      </AccordionWrapper>
    );
  };

  const renderNumSeqPerLocationBarPlot = () => {
    return (
      <AccordionWrapper
        title="Number of Sequences Per Location"
        defaultCollapsed={false}
        maxHeight={'500px'}
        helpText={
          <ul>
            <li>
              This plot shows the number of sequences{' '}
              <b>{configStore.getGroupLabel()}s</b> per location.
            </li>
            <li>
              This can be used to determine biases where large amounts of
              sequences can potentially effect how certain visualizations and
              data will appear.
            </li>
          </ul>
        }
      >
        <NumSeqPerLocationBar width={width - 200} />
      </AccordionWrapper>
    );
  };

  const renderNumSeqPerLocationLinePlot = () => {
    return (
      <AccordionWrapper
        title="Number of Sequences Per Location Over Time"
        defaultCollapsed={false}
        maxHeight={'500px'}
        helpText={
          <ul>
            <li>
              This plot shows the number of sequences{' '}
              <b>{configStore.getGroupLabel()}s</b> per location over the
              selected time frames.
            </li>
            <li>
              This can be used to determine biases where large amounts of
              sequences can potentially effect how certain visualizations and
              data will appear.
            </li>
          </ul>
        }
      >
        <NumSeqPerLocationLine width={width - 200} />
      </AccordionWrapper>
    );
  };

  return (
    <CompareGroupsTabContainer ref={ref}>
      {renderEntropyPlot()}
      {renderCooccurrencePlot()}
      {renderGroupStackPlot()}
      {renderMutationStructureViewer()}
      {renderNumSeqPerLocationBarPlot()}
      {renderNumSeqPerLocationLinePlot()}
    </CompareGroupsTabContainer>
  );
});

export default CompareGroupsTab;
