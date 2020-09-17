import React, { useEffect } from 'react';
import { ASYNC_STATES } from '../constants/UI';
import { useStores } from '../stores/connect';

const WaitForAsyncData = ({ children }) => {
  const { asyncData } = useStores();
  useEffect(() => {
    if (asyncData.status === ASYNC_STATES.UNINITIALIZED) {
      asyncData.fetchData();
    }
  });
  if (asyncData.status !== ASYNC_STATES.SUCCEEDED) {
    return <div>waiting</div>;
  }
  return children;
};

export default WaitForAsyncData;
