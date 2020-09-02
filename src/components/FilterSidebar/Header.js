import React from 'react';
import styled from 'styled-components';

import ExternalLink from '../Common/ExternalLink';
import CGLogo from '../../assets/images/cg_logo_v13.png';
import GISAIDLogo from '../../assets/images/gisaid_logo.png';

const HeaderDiv = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;
  justify-content: flex-start;
  border-bottom: 1px solid #aaa;
  flex-shrink: 0;
`;
const TitleContainer = styled.div`
  display: flex;
  flex-direction: column;
  align-items: center;
  border-bottom: 1px solid #aaa;
  padding: 5px 0px;

  background-color: #fff;

  h1 {
    font-weight: 700;
    font-size: 1.1em;
    margin: 0px;
    line-height: normal;
  }
`;

const ImageContainer = styled.div`
  margin-bottom: 2px;
  img {
    width: auto;
    height: 60px;
  }
`;

const GISAIDContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: center;

  font-weight: bold;
  font-size: 1.1em;

  margin: 5px 12px;

  a {
    display: flex;
    flex-direction: row;
    align-items: center;

    margin-left: 5px;
    img {
      height: 36px;
    }
  }
`;

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
