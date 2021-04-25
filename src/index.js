/* eslint-disable import/default */

import React from 'react';
import { render } from 'react-dom';
import * as Sentry from '@sentry/react';
import { Integrations } from '@sentry/tracing';

import App from './components/App';

require('./favicon.ico'); // Tell webpack to load favicon.ico

if (process.env.NODE_ENV === 'production' || false) {
  Sentry.init({
    dsn:
      'https://ad5d766ba9fb46d784ee1418d837605b@o561173.ingest.sentry.io/5697727',
    integrations: [new Integrations.BrowserTracing()],

    // Set tracesSampleRate to 1.0 to capture 100%
    // of transactions for performance monitoring.
    // We recommend adjusting this value in production
    tracesSampleRate: 1.0,
  });
}

render(<App />, document.getElementById('app'));
