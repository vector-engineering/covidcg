import React from 'react';

import { RouterStore, startRouter } from 'mobx-router';
import routes from '../routes';
import ObservableDataStore from './dataStore';
import ObservableUIStore from './uiStore';

export const uiStoreInstance = new ObservableUIStore();
export const dataStoreInstance = new ObservableDataStore();
export const routerInstance = new RouterStore();

export const rootStore = {
  dataStore: dataStoreInstance,
  router: routerInstance,
  uiStore: uiStoreInstance,
};

export const storesContext = React.createContext(rootStore);

export const StoreProvider = storesContext.Provider;

startRouter(routes, rootStore);
