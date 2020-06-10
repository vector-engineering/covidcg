import React from 'react';
import { storesContext } from './rootStore';

export const useStores = () => React.useContext(storesContext);

// wrap a component export with this method and the component
// will have this.props.covidStore and this.props.router
export const connect = (Component) => {
  // eslint-disable-next-line react/display-name
  return (props) => {
    return (
      <storesContext.Consumer>
        {(value) => (
          <Component
            covidStore={value.covidStore}
            router={value.router}
            uiStore={value.uiStore}
            {...props}
          />
        )}
      </storesContext.Consumer>
    );
  };
};
