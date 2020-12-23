import { appConfig } from './config';

const _metadataFields = [];
const _metadataFieldNiceNameMap = {};

Object.keys(appConfig.metadata_cols[appConfig.ingest_method]).forEach(col => {
  const colObj = appConfig.metadata_cols[appConfig.ingest_method][col];

  // Skip disabled metadata cols
  if (Object.prototype.hasOwnProperty.call(colObj, 'disabled') && colObj['disabled']) {
    return;
  }

  _metadataFields.push(col);
  _metadataFieldNiceNameMap[col] = colObj.title;
});

export const metadataFields = _metadataFields;
export const metadataFieldNiceNameMap = _metadataFieldNiceNameMap;
