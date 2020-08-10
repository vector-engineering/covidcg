import React from 'react';

import { Route } from 'mobx-router';
import HomePage from './components/Pages/HomePage';
import NotFoundPage from './components/Pages/NotFoundPage';
import { UIStoreInstance } from './stores/rootStore';

export const publicPath = '/';

const routes = {
  home: new Route({
    path: publicPath,
    component: <HomePage />,
    onEnter: (route, params, store, queryParams) => {
      if (queryParams.tab) {
        UIStoreInstance.setActiveTab(queryParams.tab);
      }
    },
  }),
  home_index: new Route({
    path: publicPath + 'index.html',
    component: <HomePage />,
    onEnter: (route, params, store, queryParams) => {
      //store.gallery.fetchImages();
      console.log(route, params, store, queryParams);
      console.log('current query params are -> ', queryParams);
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
