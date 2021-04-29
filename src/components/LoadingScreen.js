import React, { useState, useEffect } from 'react';
import { observer } from 'mobx-react';
import { Transition } from 'react-transition-group';
import { asyncDataStoreInstance } from './App';

import { ASYNC_STATES } from '../constants/defs.json';

import ExternalLink from './Common/ExternalLink';
import LoadingSpinner from './Common/LoadingSpinner';
import {
  LoadingScreenContainer,
  LogoContainer,
  ProgressContainer,
  ProgressText,
} from './LoadingScreen.styles';

import CGLogo from '../assets/images/cg_logo_v13.svg';

const duration = 1000;

const defaultStyle = {
  transition: `all ${duration}ms ease-in-out`,
  opacity: 0,
  transform: 'scale(.8, .8)',
};

const transitionStyles = {
  entering: { opacity: 0, transform: 'scale(.8, .8)' },
  entered: { opacity: 1, transform: 'scale(1, 1)' },
  exiting: { opacity: 0 },
  exited: { opacity: 0 },
};

const loadingDefaultStyles = {
  transition: `all 1000ms ease-in-out`,
  opacity: 0,
  transform: 'scale(.8, .8)',
};

const loadingTransitionStyles = {
  entering: { opacity: 0, transform: 'scale(.8, .8)' },
  entered: { opacity: 1, transform: 'scale(1, 1)' },
  exiting: { opacity: 0 },
  exited: { opacity: 0 },
};

const LoadingScreen = observer(() => {
  const [progressText, setProgressText] = useState('');

  useEffect(() => {
    if (asyncDataStoreInstance.status === ASYNC_STATES.FAILED) {
      setProgressText(
        <span>
          Error fetching data. Please refresh this page. If this error persists,
          please contact us at{' '}
          <ExternalLink href="mailto:covidcg@broadinstitute.org">
            covidcg@broadinstitute.org
          </ExternalLink>
        </span>
      );
    } else if (asyncDataStoreInstance.status !== ASYNC_STATES.SUCCEEDED) {
      setProgressText(<span>Initializing...</span>);
    }
  }, [asyncDataStoreInstance.status]);

  useEffect(() => {
    // Display message to try a refresh if this is taking too long
    // Timer set at 5 seconds
    setTimeout(() => {
      if (asyncDataStoreInstance.status === ASYNC_STATES.STARTED) {
        setProgressText(
          <span>
            Initializing...
            <br />
            If this is taking too long, try refreshing the page.
            <br />
            If the error persists, please contact us at{' '}
            <ExternalLink href="mailto:covidcg@broadinstitute.org">
              covidcg@broadinstitute.org
            </ExternalLink>
          </span>
        );
      }
    }, 10000);
  }, []);

  return (
    <LoadingScreenContainer>
      <Transition in appear enter timeout={duration}>
        {(state) => (
          <LogoContainer
            style={{
              ...defaultStyle,
              ...transitionStyles[state],
            }}
          >
            <img src={CGLogo}></img>
          </LogoContainer>
        )}
      </Transition>
      <Transition in appear enter timeout={300}>
        {(state) => (
          <ProgressContainer
            style={{
              ...loadingDefaultStyles,
              ...loadingTransitionStyles[state],
            }}
          >
            <LoadingSpinner
              visible={asyncDataStoreInstance.status !== ASYNC_STATES.FAILED}
            />
            <ProgressText>{progressText}</ProgressText>
          </ProgressContainer>
        )}
      </Transition>
    </LoadingScreenContainer>
  );
});

export default LoadingScreen;
