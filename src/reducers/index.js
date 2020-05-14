import { combineReducers } from 'redux';
import covidReducer from './covidReducer';
import { connectRouter } from 'connected-react-router'

const rootReducer = history => combineReducers({
  covid: covidReducer,
  router: connectRouter(history)
});

export default rootReducer;
