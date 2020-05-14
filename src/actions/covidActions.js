import { 
    SELECT_GENE,
    SELECT_LOCATIONS,
    SELECT_DATE_RANGE
} from '../constants/actionTypes';

// example of a thunk using the redux-thunk middleware
// export function saveFuelSavings(settings) {
//   return function (dispatch) {
//     // thunks allow for pre-processing actions, calling apis, and dispatching multiple actions
//     // in this case at this point we could call a service that would persist the fuel savings
//     return dispatch({
//       type: types.SAVE_FUEL_SAVINGS,
//       dateModified: getFormattedDateTime(),
//       settings
//     });
//   };
// }

export function selectGene(value) {
  return {
    type: SELECT_GENE,
    value: value
  };
}

export function selectLocations(selectedNodes) {
  return {
    type: SELECT_LOCATIONS,
    selectedNodes: selectedNodes
  };
}

export function selectDateRange(dateRange) {
  return {
    type: SELECT_DATE_RANGE,
    dateRange: dateRange
  }
}
