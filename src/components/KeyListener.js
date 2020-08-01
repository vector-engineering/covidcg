import React, { useEffect } from 'react';
import { useStores } from '../stores/connect';

const KeyListener = () => {
  const { UIStore } = useStores();

  const onKeyDown = (e) => {
    UIStore.setKeyPressed(e.keyCode, true);
  };

  const onKeyUp = (e) => {
    UIStore.setKeyPressed(e.keyCode, false);
  };

  useEffect(() => {
    document.addEventListener('keydown', onKeyDown, false);
    document.addEventListener('keyup', onKeyUp, false);

    return () => {
      document.removeEventListener('keydown', onKeyDown, false);
      document.removeEventListener('keyup', onKeyUp, false);
    };
  });

  return <div />;
};

export default KeyListener;
