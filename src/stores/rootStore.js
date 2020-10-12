import React from 'react';
import { RouterStore, startRouter } from 'mobx-router';
import routes from '../routes';
import ObservableDataStore from './dataStore';
import ObservableUIStore from './UIStore';
import ObservableConfigStore from './configStore';
import ObservablePlotSettingsStore from './plotSettingsStore';
import { LineageDataStore } from './lineageData';
import { SnpDataStore } from './snpData';
import { LocationDataStore } from './locationData';

class RootStore {
  router;
  UIStore;
  plotSettingsStore;

  locationDataStore;
  snpDataStore;
  lineageDataStore;

  configStore;
  dataStore;

  init() {
    this.UIStore = new ObservableUIStore();
    this.router = new RouterStore();

    this.plotSettingsStore = new ObservablePlotSettingsStore();
    this.locationDataStore = new LocationDataStore();
    this.snpDataStore = new SnpDataStore();
    this.lineageDataStore = new LineageDataStore();

    this.configStore = new ObservableConfigStore();
    this.dataStore = new ObservableDataStore();

    startRouter(routes, this);
  }
}

export const rootStoreInstance = new RootStore();

export const storesContext = React.createContext(rootStoreInstance);

export const StoreProvider = storesContext.Provider;
