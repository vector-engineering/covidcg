import { observer } from 'mobx-react';
import React, { useEffect, useRef } from 'react';
import { ASYNC_STATES } from '../constants/UI';
import { rootStoreInstance } from '../stores/rootStore';
import { asyncDataStoreInstance } from './App';

const WaitForAsyncWrapper = observer(({ children }) => {
  const renderedOnce = useRef(false);
  const renderedOnceWithAsyncData = useRef(false);

  useEffect(() => {
    if (!renderedOnce.current) {
      asyncDataStoreInstance.fetchData();
      renderedOnce.current = true;
    }
  });

  if (asyncDataStoreInstance.status !== ASYNC_STATES.SUCCEEDED) {
    return <div>waiting</div>;
  }

  if (!renderedOnceWithAsyncData.current) {
    rootStoreInstance.init();
    renderedOnceWithAsyncData.current = true;
  }

  return children;
});

export default WaitForAsyncWrapper;
