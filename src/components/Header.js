import React from 'react';
import styled from 'styled-components';
import { Link } from 'mobx-router';

import routes from '../routes';
import { useStores } from '../stores/connect';
import { version, dataDate } from '../utils/version';

const HeaderDiv = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-start;

  border-bottom: 1px solid #aaa;
`;
const TitleContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  width: 425px;

  padding-left: 30px;

  h1 {
    font-weight: 700;
    font-size: 2em;
    margin: 0px;
  }
`;
const NavLinks = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  margin-right: 30px;

  a {
    margin-right: 15px;
  }
`;

const VersionDiv = styled.div`
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  justify-content: center;

  height: 100%;
  padding-left: 10px;
  border-left: 1px solid #aaa;

  font-weight: normal;
  font-size: 0.75em;
  color: #888888;

  span.version {
    span.version-num {
      font-weight: bold;
    }
  }
  span.data-date {
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
        <span className="version">
          Version: <span className="version-num">{version}</span>
        </span>
        <span className="data-date">
          Sequences Analyzed: Up to <span className="date">{dataDate}</span>
        </span>
      </VersionDiv>
    </HeaderDiv>
  );
};

Header.displayName = 'Header';

export default Header;
