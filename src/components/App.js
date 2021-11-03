/* eslint-disable import/no-named-as-default */
import React from 'react';
import PropTypes from 'prop-types';
import { createGlobalStyle } from 'styled-components';

import { MobxRouter } from 'mobx-router';
import { rootStoreInstance, StoreProvider } from '../stores/rootStore';
import WaitForAsyncWrapper from './WaitForAsyncWrapper';
import ObservableAsyncDataStore from '../stores/asyncDataStore';
import InitialValueStore from '../stores/initialValueStore';

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
export const initialValueStoreInstance = new InitialValueStore();

const App = () => {
  return (
    <WaitForAsyncWrapper>
      <StoreProvider value={rootStoreInstance}>
        <MobxRouter store={rootStoreInstance} />
        <GlobalStyle />
      </StoreProvider>
    </WaitForAsyncWrapper>
  );
};

App.propTypes = {
  children: PropTypes.element,
};

export default App;
