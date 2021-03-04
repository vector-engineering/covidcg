import React from 'react';

import ExternalLink from '../Common/ExternalLink';
import CGLogo from '../../assets/images/cg_logo_v13.png';
import GISAIDLogo from '../../assets/images/gisaid_logo.png';

import {
  HeaderDiv,
  TitleContainer,
  ImageContainer,
  GISAIDContainer,
} from './Header.styles';

const Header = () => {
  return (
    <HeaderDiv>
      <TitleContainer>
        <ImageContainer>
          <img src={CGLogo}></img>
        </ImageContainer>
        <h1>COVID-19 CoV Genetics</h1>
      </TitleContainer>
      <GISAIDContainer>
        Enabled by data from&nbsp;
        <ExternalLink href="https://www.gisaid.org/" showIcon={false}>
          <img src={GISAIDLogo}></img>
        </ExternalLink>
      </GISAIDContainer>
    </HeaderDiv>
  );
};

Header.displayName = 'Header';

export default Header;
