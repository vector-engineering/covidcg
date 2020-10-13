import React from 'react';
import { observer } from 'mobx-react';
import { Transition } from 'react-transition-group';
import { useStores } from '../stores/connect';
import { asyncDataStoreInstance } from './App';

import { ASYNC_STATES } from '../constants/UI';

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
  const { UIStore } = useStores();

  let progressText;
  if (asyncDataStoreInstance.status !== ASYNC_STATES.SUCCEEDED) {
    progressText = 'Downloading data...';
  } else if (UIStore.caseDataState !== ASYNC_STATES.SUCCEEDED) {
    progressText = 'Processing data...';
  }

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
            <LoadingSpinner />
            <ProgressText>{progressText}</ProgressText>
          </ProgressContainer>
        )}
      </Transition>
    </LoadingScreenContainer>
  );
});

export default LoadingScreen;
