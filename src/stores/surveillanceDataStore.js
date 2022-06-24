import { action } from 'mobx';

import { asyncDataStoreInstance } from '../components/App';
import { rootStoreInstance } from './rootStore';
import { downloadBlobURL } from '../utils/download';

export class SurveillanceDataStore {
  surv_group_counts;
  surv_group_regression;

  init() {
    this.surv_group_counts = asyncDataStoreInstance.data.surv_group_counts;

    this.surv_group_regression =
      asyncDataStoreInstance.data.surv_group_regression;
  }

  @action
  downloadGroupCounts = () => {
    try {
      const blob = new Blob([JSON.stringify(this.surv_group_counts)]);
      const url = URL.createObjectURL(blob);
      downloadBlobURL(url, 'group_counts2.csv');
      rootStoreInstance.UIStore.onDownloadFinished();
    } catch (err) {
      let prefix = 'Error downloading surveillance group counts';
      if (!(typeof err.text === 'function')) {
        console.error(prefix, err);
      } else {
        err.text().then((errMsg) => {
          console.error(prefix, errMsg);
        });
      }
      rootStoreInstance.UIStore.onDownloadErr();
    }
  };
}
