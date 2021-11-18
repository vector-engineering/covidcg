import { asyncDataStoreInstance } from '../components/App';

export class SurveillanceDataStore {
  surv_group_counts;
  surv_group_regression;

  init() {
    this.surv_group_counts = asyncDataStoreInstance.data.surv_group_counts;
    this.surv_group_regression =
      asyncDataStoreInstance.data.surv_group_regression;
  }
}
