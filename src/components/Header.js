import React from 'react';
import styled from 'styled-components';
import { Link } from 'mobx-router';

import routes from '../routes';
import { useStores } from '../stores/connect';

const HeaderDiv = styled.div`
  grid-column: col2 / col3;
  grid-row: row1 / row2;

  display: grid;
  grid-template-columns: [col1] 250px [col2] auto [col3];
  grid-template-rows: [row1] auto [row2];

  border-bottom: 1px solid #aaa;
`;
const TitleContainer = styled.div`
  grid-column: col1 / col2;
  grid-row: row1 / row2;

  display: flex;
  flex-direction: row;
  align-items: center;

  padding-left: 30px;

  h1 {
    font-weight: 700;
    font-size: 2em;
    margin: 0px;
  }
`;
const NavLinks = styled.div`
  grid-column: col2 / col3;
  grid-row: row1 / row2;

  display: flex;
  flex-direction: row;
  align-items: center;

  a {
    margin-right: 15px;
  }
`;

const Header = () => {
  const { router } = useStores();
  return (
    <HeaderDiv>
      <TitleContainer>
        <h1>COVID-UI</h1>
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
    </HeaderDiv>
  );
};

Header.displayName = 'Header';

export default Header;
