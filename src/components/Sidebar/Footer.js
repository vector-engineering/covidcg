import React from 'react';
import { useStores } from '../../stores/connect';
import { version } from '../../utils/version';

import ExternalLink from '../Common/ExternalLink';

import { config } from '../../config';

import { FooterContainer, DAA, Version, SequenceMeta } from './Footer.styles';

const Footer = () => {
  const { dataStore } = useStores();

  return (
    <FooterContainer>
      {config['show_logos']['GISAID'] && (
        <DAA>
          GISAID data provided on this website are subject to GISAID’s{' '}
          <ExternalLink href="https://www.gisaid.org/registration/terms-of-use/">
            Terms and Conditions
          </ExternalLink>
        </DAA>
      )}
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
