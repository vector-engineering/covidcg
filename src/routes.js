import React from 'react';

import { Route } from 'mobx-router';
import HomePage from './components/HomePage';
import AboutPage from './components/AboutPage';
import NotFoundPage from './components/NotFoundPage';

export const publicPath = '/';

const routes = {
  home: new Route({
    path: publicPath,
    component: <HomePage />,
  }),
  home_index: new Route({
    path: publicPath + 'index.html',
    component: <HomePage />,
  }),
  about: new Route({
    path: publicPath + 'about',
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
  notFound: new Route({
    path: '*',
    component: <NotFoundPage />,
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
