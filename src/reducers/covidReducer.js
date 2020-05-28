import {
    SELECT_GENE,
    SELECT_LOCATIONS,
    SELECT_DATE_RANGE
} from '../constants/actionTypes';
import objectAssign from 'object-assign';
import initialState from './initialState';

import {
  getGene
} from '../utils/gene';

import {
  getLocationIds
} from '../utils/location';

import {
  processCaseData,
  aggCaseDataByLineage
} from '../utils/caseData';

import {
  getLineagesFromGene,
} from '../utils/lineageData';

// IMPORTANT: Note that with Redux, state should NEVER be changed.
// State is considered immutable. Instead,
// create a copy of the state passed and set new values on the copy.
// Note that I'm using Object.assign to create a copy of current state
// and update values on the copy.
export default function covidReducer(state = initialState.covid, action) {
  let newState = objectAssign({}, state);

  switch (action.type) {

    case SELECT_GENE:
      console.log('SELECT_GENE', action);

      newState.selectedGene = action.value;
      let selectedGene = getGene(action.value);

      newState.startPos = selectedGene.start;
      newState.endPos = selectedGene.end;

      // Get matching clade_ids
      let lineages = getLineagesFromGene(selectedGene);
      newState.selectedLineages = lineages;

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
      newState.caseData = processCaseData(newState.initialCaseData, newState.selectedLocationIds);

    case SELECT_DATE_RANGE:

      let { caseDataAggLineageList, changingPositions } = aggCaseDataByLineage(
        newState.caseData, newState.initialLineageData, 
        newState.startPos, newState.endPos, newState.dateRange
      );

      newState.caseDataAggLineageList = caseDataAggLineageList;
      newState.changingPositions = changingPositions;

      break;

    default:
      break;
  }

  // Recalculate case data
  return newState;
}

