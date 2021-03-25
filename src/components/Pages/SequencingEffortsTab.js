import React from 'react';
import { useStores } from '../../stores/connect';
import { observer } from 'mobx-react';
import useDimensions from 'react-use-dimensions';

import { ASYNC_STATES } from '../../constants/defs.json';

import ExternalLink from '../Common/ExternalLink';
import SequencingMapPlot from '../Vega/SequencingMapPlot';
import SkeletonElement from '../Common/SkeletonElement';
import WarningBox from '../Common/WarningBox';

import { Container, Header, Title } from './SequencingEffortsTab.styles';

const SequencingEffortsTab = observer(() => {
  const { UIStore } = useStores();
  const [ref, { width }] = useDimensions();

  if (UIStore.globalSequencingDataState === ASYNC_STATES.FAILED) {
    return (
      <>
        <div style={{ height: '20px' }} />
        <WarningBox title={'Failed to load data'} showDismissButton={false}>
          Failed to load sequencing data. Please try again by refreshing the
          page. If this error persists, please contact us at{' '}
          <ExternalLink href="mailto:covidcg@broadinstitute.org">
            covidcg@broadinstitute.org
          </ExternalLink>
        </WarningBox>
      </>
    );
  } else if (
    UIStore.globalSequencingDataState === ASYNC_STATES.UNINITIALIZED ||
    UIStore.globalSequencingDataState === ASYNC_STATES.STARTED
  ) {
    return (
      <div
        style={{
          padding: '24px',
        }}
      >
        <div
          style={{
            display: 'flex',
            flexDirection: 'row',
            alignItems: 'center',
            height: '300px',
          }}
        >
          <SkeletonElement delay={2} width={'50%'} height={300} />
          <div style={{ width: '24px' }} />
          <SkeletonElement delay={2} width={'50%'} height={300} />
        </div>
        <div
          style={{
            display: 'flex',
            flexDirection: 'row',
            alignItems: 'center',
            height: '500px',
            marginTop: '24px',
          }}
        >
          <SkeletonElement delay={2} width={'50%'} height={500} />
          <div style={{ width: '24px' }} />
          <SkeletonElement delay={2} width={'50%'} height={500} />
        </div>
      </div>
    );
  }

  return (
    <Container ref={ref}>
      <Header>
        <Title>Global sequencing coverage</Title>
        <p>
          The number of genomic sequence and associate data shared via the
          GISAID Initiative (
          <ExternalLink href="https://doi.org/10.1002/gch2.1018">
            Elbe et al, 2017, <i>Wiley Global Challenges</i>
          </ExternalLink>
          ) and case data is obtained from{' '}
          <ExternalLink href="https://github.com/CSSEGISandData/COVID-19">
            JHU CSSE COVID-19 Data
          </ExternalLink>{' '}
          (
          <ExternalLink href="https://doi.org/10.1016/S1473-3099(20)30120-1">
            Dong et al, 2020, <i>Lancet Inf Dis.</i>
          </ExternalLink>
          ). Regions with &lt;100 confirmed cases are excluded from the bar
          graphs below. Regions with &gt;20 sequences per 1000 cases are colored
          the same in the left map.
        </p>
        <p>
          <i>
            *Please note that we are currently working on cleaning the metadata
            delineating territories and countries
          </i>
        </p>
      </Header>

      <SequencingMapPlot width={width - 250} />
    </Container>
  );
});
export default SequencingEffortsTab;
