/* eslint-disable import/no-named-as-default */
import React from 'react';
import PropTypes from 'prop-types';
import { hot } from 'react-hot-loader';

import { MobxRouter } from 'mobx-router';
import { StoreProvider, rootStore } from '../stores/rootStore';

const App = () => {
  return (
    <StoreProvider value={rootStore}>
      <MobxRouter store={rootStore} />
    </StoreProvider>
  );
};

App.propTypes = {
  children: PropTypes.element,
};

export default hot(module)(App);
