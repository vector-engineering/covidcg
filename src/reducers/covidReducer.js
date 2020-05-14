import {
    SELECT_GENE,
    SELECT_LOCATIONS,
    SELECT_DATE_RANGE
} from '../constants/actionTypes';
import objectAssign from 'object-assign';
import initialState from './initialState';
import _ from 'underscore';

import {
  getGene
} from '../utils/gene';

import {
  getLocationIds
} from '../utils/location';

import {
  processCaseData,
  aggCaseDataByClade
} from '../utils/caseData';

import {
  getCladesFromGene,
} from '../utils/cladeData';

// IMPORTANT: Note that with Redux, state should NEVER be changed.
// State is considered immutable. Instead,
// create a copy of the state passed and set new values on the copy.
// Note that I'm using Object.assign to create a copy of current state
// and update values on the copy.
export default function covidReducer(state = initialState.covid, action) {
  let newState = objectAssign({}, state);

  switch (action.type) {
    // case SAVE_FUEL_SAVINGS:
    //   // For this example, just simulating a save by changing date modified.
    //   // In a real app using Redux, you might use redux-thunk and handle the async call in fuelSavingsActions.js
    //   return objectAssign({}, state, {dateModified: action.dateModified});

    // case CALCULATE_FUEL_SAVINGS:
    //   newState = objectAssign({}, state);
    //   newState[action.fieldName] = action.value;
    //   newState.necessaryDataIsProvidedToCalculateSavings = necessaryDataIsProvidedToCalculateSavings(newState);
    //   newState.dateModified = action.dateModified;

    //   if (newState.necessaryDataIsProvidedToCalculateSavings) {
    //     newState.savings = calculateSavings(newState);
    //   }

    //   return newState;

    case SELECT_GENE:
      console.log('SELECT_GENE', action);

      newState.selectedGene = action.value;
      let selectedGene = getGene(action.value);

      newState.startPos = selectedGene.start;
      newState.endPos = selectedGene.end;

      // Get matching clade_ids
      let clades = getCladesFromGene(selectedGene);
      newState.selectedClades = clades;
      //let clade_ids = _.map(clades, clade => clade.index);
      //newState.selectedCladeIds = clade_ids;
      //console.log('Clade IDs:', clade_ids);

      break;

    case SELECT_LOCATIONS:
      console.log('SELECT_LOCATIONS', action);
      newState = objectAssign({}, state);

      newState.selectedLocationIds = getLocationIds(action.selectedNodes);

      console.log('Location IDs:', newState.selectedLocationIds);
      break;

    case SELECT_DATE_RANGE:
      console.log('SELECT_DATE_RANGE', action);
      newState.dateRange = action.dateRange;

      break;

    default:
      break;
  }

  // Depending on the action, we'll have to recalculate the case data
  switch (action.type) {
    case SELECT_GENE:
    case SELECT_LOCATIONS:
      newState.caseData = processCaseData(newState.initialCaseData, newState.initialCladeData, newState.selectedLocationIds);

    case SELECT_DATE_RANGE:

      let { caseDataAggCladeList, changingPositions } = aggCaseDataByClade(
        newState.caseData, newState.initialCladeData, 
        newState.startPos, newState.endPos, newState.dateRange
      );

      newState.caseDataAggCladeList = caseDataAggCladeList;
      newState.changingPositions = changingPositions;

      break;

    default:
      break;
  }

  // Recalculate case data
  return newState;
}

