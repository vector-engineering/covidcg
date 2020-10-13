import React from 'react';
import { observer } from 'mobx-react';
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
      <LogoContainer>
        <img src={CGLogo}></img>
      </LogoContainer>
      <ProgressContainer>
        <LoadingSpinner />
        <ProgressText>{progressText}</ProgressText>
      </ProgressContainer>
    </LoadingScreenContainer>
  );
});

export default LoadingScreen;
