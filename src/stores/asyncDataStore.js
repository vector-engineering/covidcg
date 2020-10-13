import { observable, action, runInAction } from 'mobx';
import { ASYNC_STATES } from '../constants/UI';
import { Zlib } from 'zlibjs/bin/gunzip.min.js';

class ObservableAsyncDataStore {
  @observable status = ASYNC_STATES.UNINITIALIZED;
  data = {};

  @observable globalGroupCounts = [];

  @action
  async fetchData() {
    this.status = ASYNC_STATES.STARTED;
    console.log('fetch data');
    try {
      const res = await fetch(
        'https://storage.googleapis.com/ve-public/data_package.json.gz',
        {
          headers: {
            'Accept-Encoding': 'gzip',
          },
        }
      );
      // Don't try to decode data, get it in memory as bytes
      const blob = await res.blob();
      const compressed = await blob.arrayBuffer();
      // Feed into gunzip to decompress
      const gunzip = new Zlib.Gunzip(new Uint8Array(compressed));
      // Bytes to text, then text to JSON
      const data = JSON.parse(
        new TextDecoder('utf-8').decode(gunzip.decompress())
      );
      runInAction(() => {
        this.data = data;
        console.log(data);
        this.status = ASYNC_STATES.SUCCEEDED;
      });
    } catch (e) {
      console.error(e);
      runInAction(() => {
        this.status = ASYNC_STATES.FAILED;
      });
    }
  }
}

export default ObservableAsyncDataStore;
