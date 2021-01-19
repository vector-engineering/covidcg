import React, { useEffect, useState } from 'react';
// import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import styled from 'styled-components';
import _ from 'underscore';

import Button from '../Buttons/Button';
import MultiSelect from 'react-multi-select-component';

import {
  metadataFields,
  metadataFieldNiceNameMap,
} from '../../constants/metadata';

const formWidth = '160px';

const MetaFieldSelectContainer = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  padding-left: 15px;
  padding-bottom: 10px;
`;

const SelectList = styled.div`
  padding-top: 3px;
  padding-left: 10px;
`;
const SelectContainer = styled.form`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-start;

  margin-top: 3px;
  padding-right: 15px;

  label {
    margin-right: 5px;
    font-weight: normal;
  }

  .spacer {
    flex-grow: 1;
  }

  .metadata-multi-select {
    // flex-grow: 1;
    width: ${formWidth};

    .dropdown-heading {
      height: 20px;
      padding: 0 8px;

      .dropdown-heading-value {
        font-weight: normal;
        .gray {
          color: #666;
        }
      }
    }

    .dropdown-content {
      width: 250px;
      right: 0;
      padding-top: 0px;

      z-index: 2;
      .panel-content {
        border-radius: 0px;
        .select-item {
          padding: 3px 8px;
          font-weight: normal;
        }
      }
    }
  }
`;

// const PatientAgeContainer = styled.div`
//   display: flex;
//   flex-direction: row;
//   align-items: center;

//   width: ${formWidth};
//   // padding-left: 5px;

//   span {
//     font-weight: normal;
//   }
//   span.from-label {
//   }
//   span.to-label {
//     margin-left: 5px;
//   }
//   input {
//     margin-left: 5px;
//   }
// `;

const UpdateSelectionButton = styled(Button)`
  display: ${(props) => (props.show ? 'block' : 'none')};
  font-size: 1em;
  margin-left: 10px;
  margin-right: 15px;
  margin-top: 5px;
`;
UpdateSelectionButton.defaultProps = {
  show: false,
};

// Initialize field options and selections
const initialFieldOptions = {};
const initialFieldSelected = {};
metadataFields.forEach((field) => {
  initialFieldOptions[field] = [];
  initialFieldSelected[field] = [];
});

const MetaFieldSelect = observer(() => {
  const { dataStore, configStore, metadataStore } = useStores();

  const [state, setState] = useState({
    fieldOptions: initialFieldOptions,
    fieldSelected: initialFieldSelected,
    ageRange: ['', ''],
    changed: false,
  });

  const checkChanged = (fieldSelected, ageRange) => {
    // Is the current selection different than the selection in the store?
    let changed = false;
    for (let i = 0; i < Object.keys(fieldSelected).length; i++) {
      const field = Object.keys(fieldSelected)[i];

      // If this field isn't in the store's selected object, then return true
      if (
        !Object.prototype.hasOwnProperty.call(
          configStore.selectedMetadataFields,
          field
        )
      ) {
        if (fieldSelected[field].length > 0) {
          changed = true;
          break;
        } else {
          // Not defined in either, move on
          continue;
        }
      }

      // Sort both arrays so we can compare in one go
      const storeSelectedOptions = configStore.selectedMetadataFields[
        field
      ].sort((a, b) => a - b);
      // The local versions will be strings, so convert to integer IDs
      const localSelectedOptions = _.map(
        _.pluck(fieldSelected[field], 'value'),
        (a) => parseInt(a)
      );

      // If the lengths aren't the same, then break
      if (storeSelectedOptions.length !== localSelectedOptions.length) {
        changed = true;
        break;
      }

      // Check each option
      for (let j = 0; j < storeSelectedOptions.length; j++) {
        if (storeSelectedOptions[j] !== localSelectedOptions[j]) {
          changed = true;
          break;
        }
      }
      if (changed) {
        break;
      }
    }

    // Now check for age range
    // Empty values will be null in the store, and an empty string locally
    let storeAgeRange = configStore.ageRange;
    // console.log(storeAgeRange, ageRange);
    if (
      !(
        (storeAgeRange[0] === null && ageRange[0] === '') ||
        storeAgeRange[0] === parseInt(ageRange[0])
      )
    ) {
      changed = true;
    }
    if (
      !(
        (storeAgeRange[1] === null && ageRange[1] === '') ||
        storeAgeRange[1] === parseInt(ageRange[1])
      )
    ) {
      changed = true;
    }

    return changed;
  };

  // Update options when new data comes in
  useEffect(() => {
    const fieldOptions = {};
    metadataFields.forEach((field) => {
      fieldOptions[field] = [];

      if (
        !Object.prototype.hasOwnProperty.call(dataStore.metadataCounts, field)
      ) {
        return;
      }

      Object.keys(dataStore.metadataCounts[field]).forEach((option) => {
        let initialCount = dataStore.metadataCounts[field][option];
        let afterFilteringCount =
          dataStore.metadataCountsAfterFiltering[field][option];
        if (afterFilteringCount === undefined) {
          afterFilteringCount = 0;
        }
        fieldOptions[field].push({
          label:
            metadataStore.getMetadataValueFromId(field, option) +
            ' [' +
            initialCount.toString() +
            ' > ' +
            afterFilteringCount.toString() +
            ']',
          value: option,
        });
      });
    });

    setState({ ...state, fieldOptions });
  }, [dataStore.metadataCounts, dataStore.metadataCountsAfterFiltering]);

  const setSelected = (field, options) => {
    // console.log(field, options);
    let fieldSelected = Object.assign({}, state.fieldSelected);
    fieldSelected[field] = options;

    setState({
      ...state,
      fieldSelected,
      changed: checkChanged(fieldSelected, state.ageRange),
    });
  };

  // const onChangeAgeRange = (ind, e) => {
  //   // console.log(ind, e.target.value);
  //   let ageRange = state.ageRange;
  //   ageRange[ind] = e.target.value;

  //   setState({
  //     ...state,
  //     ageRange,
  //     changed: checkChanged(state.fieldSelected, ageRange),
  //   });
  // };

  // Flush changes to the store
  const updateSelectedMetadataFields = () => {
    // Only push the values, not the labels, of the selected options
    let selectedFields = Object.assign({}, state.fieldSelected);
    Object.keys(selectedFields).forEach((field) => {
      selectedFields[field] = selectedFields[field].map((option) =>
        parseInt(option.value)
      );
    });

    let ageRange = state.ageRange.slice();
    // Convert empty values to nulls, and strings to ints, for the store
    if (ageRange[0] === '') {
      ageRange[0] = null;
    } else {
      ageRange[0] = parseInt(ageRange[0]);
    }
    if (ageRange[1] === '') {
      ageRange[1] = null;
    } else {
      ageRange[1] = parseInt(ageRange[1]);
    }

    // console.log(selectedFields, ageRange);
    configStore.updateSelectedMetadataFields(selectedFields, ageRange);

    setState({
      ...state,
      changed: false
    });
  };

  // Build all of the select components
  const fieldSelects = [];
  metadataFields.forEach((field) => {
    const fieldNiceName = metadataFieldNiceNameMap[field];
    fieldSelects.push(
      <SelectContainer key={field + '_metadata_select'}>
        <label>{fieldNiceName}:</label>
        <div className="spacer"></div>
        <MultiSelect
          className={'metadata-multi-select'}
          options={state.fieldOptions[field]}
          value={state.fieldSelected[field]}
          onChange={setSelected.bind(this, field)}
          labelledBy={field}
          disableSearch={true}
          overrideStrings={{
            selectSomeItems: 'Select...',
            allItemsAreSelected: 'All items are selected.',
            selectAll: 'Select All',
            search: 'Search',
          }}
        />
      </SelectContainer>
    );
  });

  return (
    <MetaFieldSelectContainer>
      <SelectList>
        {/* <SelectContainer>
          <label>Patient Age:</label>
          <div className="spacer"></div>
          <PatientAgeContainer>
            <span className="from-label">From</span>
            <input
              type="number"
              step={1}
              min={0}
              max={120}
              value={state.ageRange[0]}
              onChange={onChangeAgeRange.bind(this, 0)}
            ></input>
            <span className="to-label">To</span>
            <input
              type="number"
              step={1}
              min={0}
              max={120}
              value={state.ageRange[1]}
              onChange={onChangeAgeRange.bind(this, 1)}
            ></input>
          </PatientAgeContainer>
        </SelectContainer> */}
        {fieldSelects}
      </SelectList>
      <UpdateSelectionButton
        show={state.changed}
        onClick={updateSelectedMetadataFields}
      >
        Update Metadata Filters
      </UpdateSelectionButton>
    </MetaFieldSelectContainer>
  );
});

export default MetaFieldSelect;
