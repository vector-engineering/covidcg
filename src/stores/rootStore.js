import React from 'react';

import { RouterStore, startRouter } from 'mobx-router';
import routes from '../routes';
import ObservableCovidStore from './covidStore';

const covidStore = new ObservableCovidStore();
const router = new RouterStore();

export const rootStore = {
  covidStore,
  router,
};

export const storesContext = React.createContext(rootStore);

export const StoreProvider = storesContext.Provider;

startRouter(routes, rootStore);
