import React from 'react';

import ExternalLink from '../Common/ExternalLink';
import CGLogo from '../../assets/images/cg_logo_v13.png';
import RSVCGLogo from '../../assets/images/rsv_pathmut_logo_v1@2x.png';
import GISAIDLogo from '../../assets/images/gisaid_logo.png';
import NCBILogo from '../../assets/images/ncbi_logo.svg';

import { config } from '../../config';

import {
  HeaderDiv,
  TitleContainer,
  ImageContainer,
  GISAIDContainer,
  NCBIContainer,
} from './Header.styles';

let logoImage;
let siteTitle;
if (config.virus === 'sars2') {
  logoImage = CGLogo;
  siteTitle = 'COVID-19 CoV Genetics';
} else if (config.virus === 'rsv') {
  logoImage = RSVCGLogo;
  siteTitle = 'RSV Pathogen Mutation DB';
}

const Header = () => {
  return (
    <HeaderDiv>
      <TitleContainer>
        <ImageContainer>
          <img src={logoImage}></img>
        </ImageContainer>
        <h1>{siteTitle}</h1>
      </TitleContainer>
      {config['show_logos']['GISAID'] && (
        <GISAIDContainer>
          Enabled by data from&nbsp;
          <ExternalLink href="https://www.gisaid.org/" showIcon={false}>
            <img src={GISAIDLogo}></img>
          </ExternalLink>
        </GISAIDContainer>
      )}
      {config['show_logos']['GenBank'] && (
        <NCBIContainer>
          <ExternalLink
            href="https://www.ncbi.nlm.nih.gov/genbank/"
            showIcon={false}
          >
            <img src={NCBILogo}></img>
          </ExternalLink>
          <span>Data from </span>
          <ExternalLink href="https://www.ncbi.nlm.nih.gov/genbank/">
            NIH GenBank
          </ExternalLink>
        </NCBIContainer>
      )}
    </HeaderDiv>
  );
};

Header.displayName = 'Header';

export default Header;
