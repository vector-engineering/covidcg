import React from 'react';
import { storesContext } from './rootStore';

export const connect = (Component) => {
  // eslint-disable-next-line react/display-name
  return (props) => {
    return (
      <storesContext.Consumer>
        {(value) => (
          <Component
            covidStore={value.covidStore}
            router={value.router}
            {...props}
          />
        )}
      </storesContext.Consumer>
    );
  };
};
