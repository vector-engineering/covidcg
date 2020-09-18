import React from 'react';

import { RouterStore, startRouter } from 'mobx-router';
import routes from '../routes';
import ObservableDataStore from './dataStore';
import ObservableUIStore from './UIStore';
import ObservableConfigStore from './configStore';
import ObservablePlotSettingsStore from './plotSettingsStore';

export const routerInstance = new RouterStore();
export const configStoreInstance = new ObservableConfigStore();
export const UIStoreInstance = new ObservableUIStore();
export const dataStoreInstance = new ObservableDataStore();
export const plotSettingsStoreInstance = new ObservablePlotSettingsStore();

export const rootStore = {
  router: routerInstance,
  configStore: configStoreInstance,
  UIStore: UIStoreInstance,
  dataStore: dataStoreInstance,
  plotSettingsStore: plotSettingsStoreInstance,
};

export const storesContext = React.createContext(rootStore);

export const StoreProvider = storesContext.Provider;

startRouter(routes, rootStore);
