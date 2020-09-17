/* eslint-disable import/no-named-as-default */
import React from 'react';
import PropTypes from 'prop-types';
import { hot } from 'react-hot-loader';
import { createGlobalStyle } from 'styled-components';

import { MobxRouter } from 'mobx-router';
import { StoreProvider, rootStore } from '../stores/rootStore';
import WaitForAsyncData from './WaitForAsyncData';

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

const App = () => {
  return (
    <StoreProvider value={rootStore}>
      <WaitForAsyncData>
        <MobxRouter store={rootStore} />
      </WaitForAsyncData>
      <GlobalStyle />
    </StoreProvider>
  );
};

App.propTypes = {
  children: PropTypes.element,
};

export default hot(module)(App);
