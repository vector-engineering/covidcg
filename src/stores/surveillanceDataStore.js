import { asyncDataStoreInstance } from '../components/App';
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
}
