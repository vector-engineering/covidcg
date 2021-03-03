import React from 'react';
import { useStores } from '../../stores/connect';
import styled from 'styled-components';

import { version } from '../../utils/version';

import ExternalLink from '../Common/ExternalLink';

const FooterContainer = styled.div`
  margin-top: auto;
  display: flex;
  flex-direction: column;
  align-items: stretch;
  background-color: #f8f8f8;

  padding: 5px 8px;
  border-top: 1px solid #ccc;

  font-size: 0.7rem;
  font-weight: normal;
  line-height: normal;
`;

const DAA = styled.div`
  margin-bottom: 5px;
  padding-bottom: 5px;
  border-bottom: 1px solid #ccc;
`;

const Version = styled.div`
  margin-bottom: 2px;
  span.version-num {
    font-weight: bold;
  }
  a {
    margin-left: 3px;
  }
`;

const SequenceMeta = styled.div`
  span.date {
    font-weight: bold;
  }
`;

const Footer = () => {
  const { dataStore } = useStores();

  return (
    <FooterContainer>
      <DAA>
        GISAID data provided on this website are subject to GISAID’s{' '}
        <ExternalLink href="https://www.gisaid.org/registration/terms-of-use/">
          Terms and Conditions
        </ExternalLink>
      </DAA>
      <Version>
        Version: <span className="version-num">{version}</span>{' '}
        <ExternalLink href="https://github.com/vector-engineering/covidcg/releases">
          (Changelog)
        </ExternalLink>
        <ExternalLink href="https://github.com/vector-engineering/covidcg">
          [GitHub]
        </ExternalLink>
      </Version>
      <SequenceMeta>
        <b>{dataStore.numSequences}</b> Sequences (Up to{' '}
        <span className="date">{dataStore.dataDate}</span>)
      </SequenceMeta>
    </FooterContainer>
  );
};

export default Footer;
