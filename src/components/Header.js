import React from 'react';
import styled from 'styled-components';
import { Link } from 'mobx-router';

import routes from '../routes';
import { useStores } from '../stores/connect';
import { version, dataDate } from '../utils/version';

const HeaderDiv = styled.div`
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  justify-content: flex-start;
  padding-top: 5px;
  padding-left: 12px;
  border-bottom: 1px solid #aaa;
`;
const TitleContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  // margin-bottom: 1px;
  h1 {
    font-weight: 700;
    font-size: 1.25em;
    margin: 0px;
    line-height: 30px;
  }
`;

const GISAIDContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;

  font-weight: normal;
  font-size: 0.85em;

  margin-bottom: 7px;

  a {
    display: flex;
    flex-direction: row;
    align-items: center;

    margin-left: 2px;
    img {
      height: 16px;
    }
  }
`;

const NavLinks = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  margin-right: 30px;
  margin-bottom: 10px;

  a {
    margin-right: 15px;
  }
`;

const VersionDiv = styled.div`
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  justify-content: center;
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
  }
  .data-date {
    line-height: normal;
    span.date {
      font-weight: bold;
    }
  }
`;

const Header = () => {
  const { router } = useStores();
  return (
    <HeaderDiv>
      <TitleContainer>
        <h1>COVID-19 CoV Genetics (CG)</h1>
      </TitleContainer>
      <GISAIDContainer>
        Powered by data from:&nbsp;
        <a
          href="https://www.gisaid.org/"
          target="_blank"
          rel="noopener noreferrer"
        >
          <img src="https://storage.googleapis.com/ve-public/covid_ui/assets/img/gisaid.png"></img>
        </a>
      </GISAIDContainer>
      <NavLinks>
        <Link router={router} route={routes.home}>
          Home
        </Link>
        <Link router={router} route={routes.about}>
          About
        </Link>
        <a
          href="https://github.com/vector-engineering/covid_ui"
          target="_blank"
          rel="noopener noreferrer"
        >
          View on GitHub
        </a>
      </NavLinks>
      <VersionDiv>
        <div className="version">
          Version: <span className="version-num">{version}</span>
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
