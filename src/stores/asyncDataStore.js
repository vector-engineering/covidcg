import { observable, action, runInAction } from 'mobx';
import { ASYNC_STATES } from '../constants/UI';

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
        'https://storage.googleapis.com/ve-public/data_package.json'
      );
      const data = await res.json();
      runInAction(() => {
        this.data = data;
        console.log(data);
        this.updateGlobalGroupCounts();
        this.status = ASYNC_STATES.SUCCEEDED;
      });
    } catch (e) {
      console.error(e);
      runInAction(() => {
        this.status = ASYNC_STATES.FAILED;
      });
    }
  }

  updateGlobalGroupCounts() {
    // Make a copy
    const processedGlobalGroupCounts = Object.assign(
      {},
      this.data.global_group_counts
    );

    // Replace integer IDs with SNP strings
    Object.keys(processedGlobalGroupCounts.dna_snp).forEach((snpId) => {
      processedGlobalGroupCounts.dna_snp[snpId.toString()] =
        processedGlobalGroupCounts.dna_snp[snpId];
    });

    Object.keys(processedGlobalGroupCounts.gene_aa_snp).forEach((snpId) => {
      processedGlobalGroupCounts.gene_aa_snp[snpId.toString()] =
        processedGlobalGroupCounts.gene_aa_snp[snpId];
    });

    Object.keys(processedGlobalGroupCounts.protein_aa_snp).forEach((snpId) => {
      processedGlobalGroupCounts.protein_aa_snp[snpId.toString()] =
        processedGlobalGroupCounts.protein_aa_snp[snpId];
    });

    this.globalGroupCounts = processedGlobalGroupCounts;
  }
}

export default ObservableAsyncDataStore;
