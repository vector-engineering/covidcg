/* eslint-disable import/no-named-as-default */
import React, { useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import { hot } from 'react-hot-loader';
import { createGlobalStyle } from 'styled-components';

import { MobxRouter } from 'mobx-router';
import ObservableAsyncDataStore from '../stores/asyncDataStore';
import { ASYNC_STATES } from '../constants/UI';
import { rootStoreInstance, StoreProvider } from '../stores/rootStore';

const GlobalStyle = createGlobalStyle`
  body {
    font: 14px 'Helvetica Neue', Helvetica, Arial, sans-serif;
    line-height: 1.4em;
    color: #4d4d4d;
    min-width: 230px;
    // max-width: 550px;
    //padding: 0px 20px;
    margin: 0 auto;
    font-smoothing: antialiased;
    font-weight: 500;
  }
`;

export const asyncDataStoreInstance = new ObservableAsyncDataStore();

const App = () => {
  if (asyncDataStoreInstance.status !== ASYNC_STATES.SUCCEEDED) {
    return <div>waiting</div>;
  }

  const renderedOnce = useRef(false);
  useEffect(() => {
    if (!renderedOnce.current) {
      rootStoreInstance.init();
      renderedOnce.current = true;
    }
  });
  return (
    <StoreProvider value={rootStoreInstance}>
      <MobxRouter store={rootStoreInstance} />
      <GlobalStyle />
    </StoreProvider>
  );
};

App.propTypes = {
  children: PropTypes.element,
};

export default hot(module)(App);
