import Worker from './example.worker.js';

const worker = new Worker();

worker.onmessage = (e) => {
  console.log('on message', e.data);
};

worker.postMessage('hello');
