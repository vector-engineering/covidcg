import React from 'react';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';
import useDimensions from 'react-use-dimensions';

import KBD from '../Common/KBD';
import TabIndicator from '../Common/TabIndicator';
import SelectBoxText from '../Common/SelectBoxText';
import AccordionWrapper from '../Common/AccordionWrapper';

import VegaStackedBars from '../Vega/GroupStackPlot';
import DataTableContainer from '../Table/DataTableContainer';
import LocationGroupPlot from '../Vega/LocationGroupPlot';
import EntropyPlot from '../Vega/EntropyPlot';
import CooccurrencePlot from '../Vega/CooccurrencePlot';
import AppStatusBox from '../Vega/AppStatusBox';
// import AcknowledgementsTable from '../Table/AcknowledgementsTable';

import { GROUP_SNV, DNA_OR_AA, TABS } from '../../constants/defs.json';

const GroupTabContainer = styled.div`
  padding-top: 10px;
`;

const GroupTab = observer(() => {
  const { configStore } = useStores();
  const [ref, { width }] = useDimensions();

  const renderAppStatusBox = () => {
    return (
      <AccordionWrapper
        title="Status"
        defaultCollapsed={false}
        maxHeight={'300px'}
      >
        <AppStatusBox />
      </AccordionWrapper>
    );
  };

  const renderEntropyPlot = () => {
    if (configStore.groupKey !== GROUP_SNV) {
      return null;
    }

    return (
      <AccordionWrapper
        title={`${configStore.getGroupLabel()} Frequencies`}
        defaultCollapsed={false}
        maxHeight={'400px'}
        helpText={
          <ul>
            <li>
              This plot shows the frequency of {configStore.getGroupLabel()}s
              along the WIV04 hCoV19 genome (
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
              <TabIndicator tab={TABS.TAB_LOCATION}>
                Compare Locations
              </TabIndicator>{' '}
              tab.
            </li>
          </ul>
        }
      >
        <EntropyPlot width={width - 150} />
      </AccordionWrapper>
    );
  };

  const renderCooccurrencePlot = () => {
    if (configStore.groupKey !== GROUP_SNV) {
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
              SNVs are counted per-co-occurrence, i.e., one matching sequence
              may count towards multiple SNVs. Only the relative counts between
              co-occurring SNVs within the same bar should be interpreted – the
              sum of all SNV counts per bar is not meaningful.
            </li>
            <li>
              <i>Click</i> on a SNV bar, or on a SNV y-axis label, to select or
              deselect a SNV. If the SNV is already selected, then clicking will
              deselect it. If it is not already selected, then clicking on the
              SNV will add it to the existing selection.
            </li>
            <li>
              SNV frequencies can be shown as raw counts{' '}
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
            {configStore.groupKey === GROUP_SNV && (
              <li>
                In <b>SNV Mode</b>, sequences are split into two groups: those
                that have the SNV (or combination of SNVs), and those that
                don&apos;t.
              </li>
            )}
            <li>
              <i>Click</i> to select one, or hold <KBD>Shift</KBD> and{' '}
              <i>click</i> to select multiple {configStore.getGroupLabel()}s.
            </li>
            <li>
              Selected {configStore.getGroupLabel()}s will be highlighted in the
              legend and table below, as well as in the{' '}
              <TabIndicator tab={TABS.TAB_LOCATION}>
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
        {configStore.groupKey !== GROUP_SNV && (
          <LocationGroupPlot width={width - 250} />
        )}
      </AccordionWrapper>
    );
  };

  const renderDataTable = () => {
    return (
      <AccordionWrapper
        title="Table"
        defaultCollapsed={false}
        maxHeight={'1200px'}
        helpText={
          <ul>
            <li>
              The table shows the counts and associated mutations of each{' '}
              <b>{configStore.getGroupLabel()}</b>.
            </li>
            <li>
              <i>Click</i> on a column to sort rows by that column. <i>Click</i>{' '}
              on the same column again to change the sort order.
            </li>
            <li>
              <i>Click</i> a table row to select one, or hold <KBD>Shift</KBD>{' '}
              and
              <i>click</i> to select multiple {configStore.getGroupLabel()}s.
            </li>
            <li>
              Selected {configStore.getGroupLabel()}s will be highlighted in the
              legend and plot, as well as in the{' '}
              <TabIndicator tab={TABS.TAB_LOCATION}>
                Compare Locations
              </TabIndicator>{' '}
              tab.
            </li>
            <li>
              {configStore.dnaOrAa === DNA_OR_AA.DNA ? 'Bases' : 'Residues'} can
              be colored by comparing them to the WIV04 reference sequence{' '}
              <SelectBoxText>Comparison to Reference</SelectBoxText>, by a set
              color code{' '}
              <SelectBoxText>
                {configStore.dnaOrAa === DNA_OR_AA.DNA ? '4' : '20'}-Color Code
              </SelectBoxText>
              , or other predefined color schemes.
            </li>
            <li>
              Comparisons to the reference sequence can be a{' '}
              <SelectBoxText>Match</SelectBoxText>, or a mismatch{' '}
              <SelectBoxText>Don&apos;t Match</SelectBoxText>.
            </li>
            <li>
              Color bases that satisfy the reference sequence matching condition
              with a set color <SelectBoxText>Yellow</SelectBoxText>,{' '}
              <SelectBoxText>Green</SelectBoxText>, etc, or with other
              predefined color schemes.
            </li>
          </ul>
        }
      >
        <DataTableContainer />
      </AccordionWrapper>
    );
  };

  // const renderAckTable = () => {
  //   return (
  //     <AccordionWrapper
  //       title="acknowledgements"
  //       defaultCollapsed={true}
  //       maxHeight={'1200px'}
  //     >
  //       <AcknowledgementsTable />
  //     </AccordionWrapper>
  //   );
  // };

  return (
    <GroupTabContainer ref={ref}>
      {renderAppStatusBox()}
      {renderEntropyPlot()}
      {renderCooccurrencePlot()}
      {renderGroupStackPlot()}
      {renderDataTable()}
      {/* renderAckTable() */}
    </GroupTabContainer>
  );
});

export default GroupTab;
