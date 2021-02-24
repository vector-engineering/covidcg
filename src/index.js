/* eslint-disable import/default */

import React from 'react';
import { render } from 'react-dom';
import App from './components/App';

require('./favicon.ico'); // Tell webpack to load favicon.ico

render(
  <AppContainer>
    <App />
  </AppContainer>,
  document.getElementById('app')
);
