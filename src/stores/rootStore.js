import React from 'react';

import { RouterStore, startRouter } from 'mobx-router';
import routes from '../routes';
import ObservableCovidStore from './covidStore';
import UiStore from './uiStore';

const covidStore = new ObservableCovidStore();
const router = new RouterStore();
const uiStore = new UiStore();

export const rootStore = {
  covidStore,
  router,
  uiStore,
};

export const storesContext = React.createContext(rootStore);

export const StoreProvider = storesContext.Provider;

startRouter(routes, rootStore);
