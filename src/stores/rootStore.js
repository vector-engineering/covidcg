import React from 'react';
import { RouterStore, startRouter } from 'mobx-router';
import routes from '../routes';
import ObservableDataStore from './dataStore';
import ObservableUIStore from './UIStore';
import ObservableConfigStore from './configStore';
import ObservablePlotSettingsStore from './plotSettingsStore';
import { lineageDataStoreInstance } from '../utils/lineageData';
import { snpDataStoreInstance } from '../utils/snpData';

class RootStore {
  router;
  configStore;
  UIStore;
  dataStore;
  plotSettingsStore;

  init() {
    this.UIStore = new ObservableUIStore();
    this.router = new RouterStore();
    this.configStore = new ObservableConfigStore();
    this.dataStore = new ObservableDataStore();
    this.plotSettingsStore = new ObservablePlotSettingsStore();
    startRouter(routes, this);
    lineageDataStoreInstance.init();
    snpDataStoreInstance.init();
  }
}

export const rootStoreInstance = new RootStore();

export const storesContext = React.createContext(rootStoreInstance);

export const StoreProvider = storesContext.Provider;

export const blah = '12';
