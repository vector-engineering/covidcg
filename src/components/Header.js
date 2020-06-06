import React from 'react';

import { Link } from 'mobx-router';

import routes from '../routes';
import { useStores } from '../stores/connect';

const Header = () => {
  const { router } = useStores();
  return (
    <div className="header">
      <div className="title-container">
        <h1>COVID-UI</h1>
      </div>
      <div className="nav-links">
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
      </div>
    </div>
  );
};

Header.displayName = 'Header';

export default Header;
