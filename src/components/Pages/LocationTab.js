import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';
import { observer } from 'mobx-react';

import ExternalLink from '../Common/ExternalLink';
import KBD from '../Common/KBD';
import AccordionWrapper from '../Common/AccordionWrapper';

import VegaLegend from '../Vega/VegaLegend';
import LocationGroupPlot from '../Vega/LocationGroupPlot';
import LocationDatePlot from '../Vega/LocationDatePlot';

import { GROUP_KEYS } from '../../constants/config';

const LocationTabContainer = styled.div`
  padding-top: 10px;
`;

const LocationTab = observer(({ width }) => {
  const { configStore } = useStores();

  return (
    <LocationTabContainer>
      <AccordionWrapper
        title="Legend"
        defaultCollapsed={false}
        maxHeight={'500px'}
        helpText={
          <p>
            Items in the legend represent <b>{configStore.getGroupLabel()}s</b>.
            Click to select one, or hold <KBD>Shift</KBD> and click to select
            multiple {configStore.getGroupLabel()}s. Sequence counts of the
            selected {configStore.getGroupLabel()}s will be compared between
            locations in the plot below.
            {configStore.groupKey === GROUP_KEYS.GROUP_LINEAGE && (
              <>
                {' '}
                <ExternalLink href="https://cov-lineages.org/descriptions.html">
                  (Lineage Descriptions)
                </ExternalLink>
              </>
            )}
          </p>
        }
      >
        <VegaLegend />
      </AccordionWrapper>
      <AccordionWrapper
        title="Location-Date Plot"
        defaultCollapsed={false}
        maxHeight={'800px'}
        helpText={
          <p>
            This plot compares the sequence counts or percentages, of the
            selected <b>{configStore.getGroupLabel()}s</b>, between the selected
            locations. Click to highlight one, or hold <KBD>Shift</KBD> and
            click to highlight multiple locations. Highlighted locations will be
            shown in the plot below as well.
          </p>
        }
      >
        <LocationDatePlot width={width - 200} />
      </AccordionWrapper>
      <AccordionWrapper
        title="Location-Group Plot"
        defaultCollapsed={false}
        maxHeight={'500px'}
        helpText={
          <p>
            This plot shows the cumulative proportion of{' '}
            <b>{configStore.getGroupLabel()}s</b> per location. Click to select
            one, or hold <KBD>Shift</KBD> and click to select multiple{' '}
            {configStore.getGroupLabel()}s. Sequences from the selected{' '}
            {configStore.getGroupLabel()}s will be shown in the plot above.
          </p>
        }
      >
        <LocationGroupPlot width={width - 300} />
      </AccordionWrapper>
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
