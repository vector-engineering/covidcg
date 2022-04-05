import React from 'react';
import { RouterStore, startRouter } from 'mobx-router';
import routes from '../routes';
import { DataStore } from './dataStore';
import { UIStore } from './UIStore';
import { ConfigStore } from './configStore';
import { PlotSettingsStore } from './plotSettingsStore';
import { MutationDataStore } from './mutationData';
import { LocationDataStore } from './locationDataStore';
import { MetadataStore } from './metadataStore';
import { GlobalSequencingDataStore } from './globalSequencingData';
import { GroupDataStore } from './groupDataStore';
import { ExampleStore } from './exampleStore';
import { SurveillanceDataStore } from './surveillanceDataStore';

class RootStore {
  UIStore;
  router;

  plotSettingsStore;
  metadataStore;
  locationDataStore;
  mutationDataStore;

  configStore;
  dataStore;

  globalSequencingDataStore;
  groupDataStore;

  surveillanceDataStore;
  exampleStore;

  constructor() {
    this.UIStore = new UIStore();
    this.router = new RouterStore();

    this.plotSettingsStore = new PlotSettingsStore();
    this.metadataStore = new MetadataStore();
    this.locationDataStore = new LocationDataStore();
    this.mutationDataStore = new MutationDataStore();

    this.configStore = new ConfigStore();
    this.dataStore = new DataStore();

    this.globalSequencingDataStore = new GlobalSequencingDataStore();
    this.groupDataStore = new GroupDataStore();

    this.exampleStore = new ExampleStore();
    this.surveillanceDataStore = new SurveillanceDataStore();
  }

  init() {
    // Initialize all stores
    this.UIStore.init();
    startRouter(routes, this);

    this.plotSettingsStore.init();
    this.metadataStore.init();
    this.locationDataStore.init();

    this.configStore.init();

    this.mutationDataStore.init();

    this.globalSequencingDataStore.init();
    this.groupDataStore.init();

    this.surveillanceDataStore.init();

    this.dataStore.init();

    this.exampleStore.init();

    const urlParams = new URLSearchParams(window.location.search);

    if (urlParams.get('tab')) {
      this.UIStore.setActiveTab(urlParams.get('tab'));
    }
  }
}

export const rootStoreInstance = new RootStore();
export const storesContext = React.createContext(rootStoreInstance);
export const StoreProvider = storesContext.Provider;
