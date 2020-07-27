import React from 'react';
import styled from 'styled-components';

import ExternalLink from '../Common/ExternalLink';
import CGLogo from '../../assets/images/cg_logo_v13.png';
import { version, dataDate } from '../../utils/version';

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
  align-items: flex-start;
  margin-bottom: 5px;
  border-bottom: 1px solid #aaa;

  background-color: #fff;

  h1 {
    font-weight: 700;
    font-size: 1.25em;
    margin: 0px;
    line-height: 30px;
    margin-left: 12px;
  }

  img {
    width: auto;
    margin-left: auto;
  }
`;

const GISAIDContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;

  font-weight: normal;
  font-size: 1em;

  margin-bottom: 7px;
  margin-left: 12px;

  a {
    display: flex;
    flex-direction: row;
    align-items: center;

    margin-left: 2px;
    img {
      height: 24px;
    }
  }
`;

const VersionDiv = styled.div`
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  justify-content: center;
  margin-left: 12px;
  margin-bottom: 7px;

  // height: 30px;
  // border-top: 1px solid #ccc;

  font-weight: normal;
  font-size: 0.75em;
  color: #666666;

  .version {
    line-height: normal;
    margin-bottom: 2px;
    span.version-num {
      font-weight: bold;
    }

    a {
      margin-left: 3px;
    }
  }
  .data-date {
    line-height: normal;
    span.date {
      font-weight: bold;
    }
  }
`;

const Header = () => {
  return (
    <HeaderDiv>
      <TitleContainer>
        <img height={78} src={CGLogo}></img>
        <h1>COVID-19 CoV Genetics</h1>
      </TitleContainer>
      <GISAIDContainer>
        SARS-CoV-2 sequences from:&nbsp;
        <ExternalLink href="https://www.gisaid.org/">
          <img src="https://storage.googleapis.com/ve-public/covid_ui/assets/img/gisaid.png"></img>
        </ExternalLink>
      </GISAIDContainer>
      <VersionDiv>
        <div className="version">
          Version: <span className="version-num">{version}</span>{' '}
          <ExternalLink href="https://github.com/vector-engineering/COVID19-CG/releases">
            (Changelog)
          </ExternalLink>
          <ExternalLink href="https://github.com/vector-engineering/COVID19-CG">
            [GitHub]
          </ExternalLink>
        </div>
        <div className="data-date">
          Sequences Analyzed: Up to <span className="date">{dataDate}</span>
        </div>
      </VersionDiv>
    </HeaderDiv>
  );
};

Header.displayName = 'Header';

export default Header;
