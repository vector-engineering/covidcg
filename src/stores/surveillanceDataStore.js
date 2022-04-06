import { action } from 'mobx';

import { asyncDataStoreInstance } from '../components/App';
import { rootStoreInstance } from './rootStore';
import { downloadBlobURL } from '../utils/download';
import { intToISO, ISOToInt } from '../utils/date';
import { config } from '../config';

const today = intToISO(new Date().getTime());
const lastNDays = 30 * 12 * 5;

export class SurveillanceDataStore {
  surv_group_counts;
  surv_group_regression;

  init() {
    this.surv_group_counts = asyncDataStoreInstance.data.surv_group_counts.filter(
      (elem) => {
        return config.virus === 'sars2'
          ? elem
          : elem.collection_week >=
              ISOToInt(today) - lastNDays * 24 * 60 * 60 * 1000;
      }
    );

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
