import React from 'react';

import { RouterStore, startRouter } from 'mobx-router';
import routes from '../routes';
import ObservableCovidStore from './covidStore';
import UiStore from './uiStore';

export const uiStoreInstance = new UiStore();
export const covidStoreInstance = new ObservableCovidStore();
export const routerInstance = new RouterStore();

export const rootStore = {
  covidStore: covidStoreInstance,
  router: routerInstance,
  uiStore: uiStoreInstance,
};

export const storesContext = React.createContext(rootStore);

export const StoreProvider = storesContext.Provider;

startRouter(routes, rootStore);
