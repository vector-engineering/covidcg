import React from 'react';
import { RouterStore, startRouter } from 'mobx-router';
import routes from '../routes';
import ObservableDataStore from './dataStore';
import ObservableUIStore from './UIStore';
import ObservableConfigStore from './configStore';
import ObservablePlotSettingsStore from './plotSettingsStore';

console.log('yo');

class RootStore {
  router;
  configStore;
  UIStore;
  dataStore;
  plotSettingsStore;

  constructor() {
    this.yo = 'yo';
  }
  init() {
    this.router = new RouterStore();
    this.configStore = new ObservableConfigStore();
    this.UIStore = new ObservableUIStore();
    this.dataStore = new ObservableDataStore();
    this.plotSettingsStore = new ObservablePlotSettingsStore();
    startRouter(routes, this);
  }
}

export const rootStoreInstance = new RootStore();

export const storesContext = React.createContext(rootStoreInstance);

export const StoreProvider = storesContext.Provider;

export const blah = '12';
