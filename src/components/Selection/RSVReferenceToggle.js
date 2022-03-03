import React from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import { RadioForm } from './GroupBySelect.styles';

// Used in the Lineage Reports tab for RSV
const RSVReferenceToggle = observer(() => {
  const { configStore, groupDataStore } = useStores();

  const defaultGenotypes = { A: 'ON1', B: 'BAIX' };

  const onReferenceChange = (event) => {
    configStore.applyPendingChanges({ selectedReference: event.target.value });
    groupDataStore.updateSelectedGroups([defaultGenotypes[event.target.value]]);
  };

  const renderRSVRefSelect = () => {
    const selectedReference = configStore.selectedReference;
    return (
      <div className="radio-row">
        <div className="radio-item">
          <input
            type="radio"
            id="RSVAChoice"
            name="rsvAorB"
            value="A"
            checked={selectedReference === 'A'}
            onChange={onReferenceChange}
          ></input>
          <label htmlFor="RSVAChoice">RSV-A</label>
        </div>
        <div className="radio-item">
          <input
            type="radio"
            id="RSVBChoice"
            name="rsvAorB"
            value="B"
            checked={selectedReference === 'B'}
            onChange={onReferenceChange}
          ></input>
          <label htmlFor="RSVAChoice">RSV-B</label>
        </div>
      </div>
    );
  };

  return (
    <RadioForm direction="row">
      <span className="form-title">Reference Sequence</span>
      {renderRSVRefSelect()}
    </RadioForm>
  );
});

export default RSVReferenceToggle;
