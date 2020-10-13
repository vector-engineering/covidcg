import React, { useEffect, useRef } from 'react';
import { observer } from 'mobx-react';
import { ASYNC_STATES } from '../constants/UI';
import { rootStoreInstance } from '../stores/rootStore';
import { asyncDataStoreInstance } from './App';

import LoadingScreen from './LoadingScreen';

const WaitForAsyncWrapper = observer(({ children }) => {
  const renderedOnce = useRef(false);
  const renderedOnceWithAsyncData = useRef(false);

  useEffect(() => {
    if (!renderedOnce.current) {
      asyncDataStoreInstance.fetchData();
      renderedOnce.current = true;
    }
  });

  if (
    asyncDataStoreInstance.status === ASYNC_STATES.SUCCEEDED &&
    !renderedOnceWithAsyncData.current
  ) {
    rootStoreInstance.init();
    renderedOnceWithAsyncData.current = true;
  }

  if (
    asyncDataStoreInstance.status !== ASYNC_STATES.SUCCEEDED ||
    rootStoreInstance.UIStore.caseDataState !== ASYNC_STATES.SUCCEEDED
  ) {
    return <LoadingScreen />;
  }

  return children;
});

export default WaitForAsyncWrapper;
