import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';

import ExternalLink from '../Common/ExternalLink';
import SequencingMapPlot from '../Vega/SequencingMapPlot';

const Container = styled.div`
  padding-top: 20px;
  overflow-x: hidden;
`;

const Header = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  width: 100%;

  padding: 0px 20px;

  p {
    font-weight: normal;
    line-height: normal;
    max-width: 800px;
    margin: 5px 0px;
  }
`;
const Title = styled.h2`
  font-size: 1.5em;
  font-weight: bold;

  margin: 0px;
  margin-bottom: 10px;
`;

const SequencingEffortsTab = ({ width }) => {
  return (
    <Container>
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

      <SequencingMapPlot width={width - 150} />
    </Container>
  );
};
SequencingEffortsTab.propTypes = {
  width: PropTypes.number,
};
SequencingEffortsTab.defaultProps = {
  width: 100,
};

export default SequencingEffortsTab;
