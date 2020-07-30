import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';
import { observer } from 'mobx-react';
import { asyncStates } from '../../stores/uiStore';

import KBD from '../Common/KBD';
import TabIndicator from '../Common/TabIndicator';
import AccordionTitle from '../Common/AccordionTitle';
import AccordionWrapper from '../Common/AccordionWrapper';
import SkeletonElement from '../Common/SkeletonElement';
import LoadingSpinner from '../Common/LoadingSpinner';

import VegaLegend from '../Vega/VegaLegend';
import LocationGroupPlot from '../Vega/LocationGroupPlot';
import LocationDatePlot from '../Vega/LocationDatePlot';

const LocationTabContainer = styled.div``;

const HelpText = styled.div`
  margin-bottom: 5px;

  font-weight: normal;
  font-size: 0.9em;
  color: #666;
  line-height: normal;
  p {
    margin-top: 0px;
    margin-bottom: 3px;
  }
`;

const LocationDatePlotWrapper = observer(({ width }) => {
  const { covidStore, uiStore } = useStores();
  if (uiStore.caseDataState === asyncStates.STARTED) {
    return (
      <div
        style={{
          paddingTop: '12px',
          paddingRight: '24px',
          paddingLeft: '12px',
          paddingBottom: '24px',
        }}
      >
        <SkeletonElement delay={2} height={300}>
          <LoadingSpinner />
        </SkeletonElement>
      </div>
    );
  } else {
    return (
      <AccordionWrapper
        title={<AccordionTitle>Location-Date Plot</AccordionTitle>}
        defaultCollapsed={false}
        maxHeight={'800px'}
      >
        <HelpText>
          <p>
            This plot compares the sequence counts or percentages, of the
            selected <b>{covidStore.getGroupLabel()}s</b>, between the selected
            locations. Click to highlight one, or hold <KBD>Shift</KBD> and
            click to highlight multiple locations. Highlighted locations will be
            shown in the plot below as well.
          </p>
        </HelpText>
        <LocationDatePlot width={width - 200} />
      </AccordionWrapper>
    );
  }
});
LocationDatePlotWrapper.propTypes = {
  width: PropTypes.number,
};

const LocationGroupPlotWrapper = observer(({ width }) => {
  const { covidStore, uiStore } = useStores();
  if (uiStore.caseDataState === asyncStates.STARTED) {
    return (
      <div
        style={{
          paddingTop: '12px',
          paddingRight: '24px',
          paddingLeft: '12px',
          paddingBottom: '24px',
        }}
      >
        <SkeletonElement delay={2} height={300}>
          <LoadingSpinner />
        </SkeletonElement>
      </div>
    );
  } else {
    return (
      <AccordionWrapper
        title={<AccordionTitle>Location-Group Plot</AccordionTitle>}
        defaultCollapsed={false}
        maxHeight={'500px'}
      >
        <HelpText>
          <p>
            This plot shows the cumulative proportion of{' '}
            <b>{covidStore.getGroupLabel()}s</b> per location. Click to select
            one, or hold <KBD>Shift</KBD> and click to select multiple{' '}
            {covidStore.getGroupLabel()}s. Sequences from the selected{' '}
            {covidStore.getGroupLabel()}s will be shown in the plot above.
          </p>
        </HelpText>
        <LocationGroupPlot width={width - 300} />
      </AccordionWrapper>
    );
  }
});
LocationGroupPlotWrapper.propTypes = {
  width: PropTypes.number,
};

const LocationTab = observer(({ width }) => {
  const { covidStore } = useStores();

  return (
    <LocationTabContainer>
      <AccordionWrapper
        title={<AccordionTitle>Legend</AccordionTitle>}
        defaultCollapsed={false}
        maxHeight={'500px'}
      >
        <HelpText>
          <p>
            Items in the legend represent <b>{covidStore.getGroupLabel()}s</b>.
            Click to select one, or hold <KBD>Shift</KBD> and click to select
            multiple {covidStore.getGroupLabel()}s. Sequence counts of the
            selected {covidStore.getGroupLabel()}s will be compared between
            locations in the plot below.
          </p>
        </HelpText>
        <VegaLegend />
      </AccordionWrapper>

      <LocationDatePlotWrapper width={width} />
      <LocationGroupPlotWrapper width={width} />
    </LocationTabContainer>
  );
});
LocationTab.propTypes = {
  width: PropTypes.number,
};
LocationTab.defaultProps = {
  width: 100,
};

export default LocationTab;
