import React from 'react';
import { storesContext } from './rootStore';

export const useStores = () => React.useContext(storesContext);

// wrap a component export with this method and the component
// will have this.props.dataStore and this.props.router
export const connect = (Component) => {
  // eslint-disable-next-line react/display-name
  return (props) => {
    return (
      <storesContext.Consumer>
        {(value) => (
          <Component
            dataStore={value.dataStore}
            router={value.router}
            UIStore={value.UIStore}
            configStore={value.configStore}
            plotSettingsStore={value.plotSettingsStore}
            locationDataStore={value.locationDataStore}
            lineageDataStore={value.lineageDataStore}
            snpDataStore={value.snpDataStore}
            metadataStore={value.metadataStore}
            globalSequencingDataStore={value.globalSequencingDataStore}
            {...props}
          />
        )}
      </storesContext.Consumer>
    );
  };
};
