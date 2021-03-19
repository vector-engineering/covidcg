import React, { useEffect, useState } from 'react';
// import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { ASYNC_STATES } from '../../constants/defs.json';

import MultiSelect from 'react-multi-select-component';
import QuestionButton from '../Buttons/QuestionButton';
import SkeletonElement from '../Common/SkeletonElement';

import {
  MetaFieldSelectContainer,
  SelectList,
  SelectContainer,
} from './MetaFieldSelect.styles';

import {
  metadataFields,
  metadataFieldNiceNameMap,
} from '../../constants/metadata';

// Initialize field options and selections
const initialFieldOptions = {};
const initialFieldSelected = {};
metadataFields.forEach((field) => {
  initialFieldOptions[field] = [];
  initialFieldSelected[field] = [];
});

const MetaFieldSelect = observer(
  ({ selectedMetadataFields, updateSelectedMetadataFields }) => {
    const { dataStore, metadataStore, UIStore } = useStores();

    const [state, setState] = useState({
      fieldOptions: initialFieldOptions,
      fieldSelected: initialFieldSelected,
      ageRange: ['', ''],
    });

    const updateOptions = () => {
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
          fieldOptions[field].push({
            label:
              metadataStore.getMetadataValueFromId(field, option) +
              ' [' +
              initialCount.toString() +
              ']',
            value: option,
          });
        });
      });

      setState({ ...state, fieldOptions });
    };

    // Update options when new data comes in
    useEffect(() => {
      updateOptions();
    }, [dataStore.metadataCounts]);

    const setSelected = (field, options) => {
      // console.log(field, options);

      let selectedFields = Object.assign({}, selectedMetadataFields);
      selectedFields[field] = options;

      updateSelectedMetadataFields(field, options);
    };

    // When the metadata fields finish loading, trigger and update
    useEffect(() => {
      if (UIStore.metadataFieldState === ASYNC_STATES.SUCCEEDED) {
        updateOptions();
      }
    }, [UIStore.metadataFieldState]);

    // const onChangeAgeRange = (ind, e) => {
    //   // console.log(ind, e.target.value);
    //   let ageRange = state.ageRange;
    //   ageRange[ind] = e.target.value;

    //   // Convert empty values to nulls, and strings to ints, for the store
    //   if (ageRange[0] === '') {
    //     ageRange[0] = null;
    //   } else {
    //     ageRange[0] = parseInt(ageRange[0]);
    //   }
    //   if (ageRange[1] === '') {
    //     ageRange[1] = null;
    //   } else {
    //     ageRange[1] = parseInt(ageRange[1]);
    //   }

    //   updateAgeRange(ageRange);
    // };

    if (UIStore.metadataFieldState !== ASYNC_STATES.SUCCEEDED) {
      return (
        <div
          style={{
            paddingTop: '12px',
            paddingRight: '24px',
            paddingLeft: '12px',
            paddingBottom: '24px',
          }}
        >
          <SkeletonElement delay={2} height={100}></SkeletonElement>
        </div>
      );
    }

    // Build all of the select components
    const fieldSelects = [];
    metadataFields.forEach((field) => {
      const fieldNiceName = metadataFieldNiceNameMap[field];
      fieldSelects.push(
        <SelectContainer key={field + '_metadata_select'}>
          <label>{fieldNiceName}:</label>
          <MultiSelect
            className={'metadata-multi-select'}
            options={state.fieldOptions[field]}
            value={selectedMetadataFields[field]}
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
        <span className="title">
          Metadata Filters
          <QuestionButton
            data-tip='<p>By default, no filtering is applied on sequence metadata (Default is select all)</p><p>Metadata is dependent on the data submitter, so many fields may be missing and marked as "Unknown".</p>'
            data-html={true}
            data-place="left"
            data-for="main-tooltip"
          />
        </span>
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
      </MetaFieldSelectContainer>
    );
  }
);

export default MetaFieldSelect;
