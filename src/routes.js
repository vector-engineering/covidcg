import React from 'react';

import { Route } from 'mobx-router';
import HomePage from './components/HomePage';
import AboutPage from './components/AboutPage';

const routes = {
  home: new Route({
    path: '/',
    component: <HomePage />,
  }),
  about: new Route({
    path: '/about',
    component: <AboutPage />,
    onEnter: () => {
      // we could do stuff here
    },
    beforeExit: () => {
      // or here
    },
    // eslint-disable-next-line no-unused-vars
    onParamsChange: (route, params, store) => {
      //or here
    },
  }),
};
export default routes;
