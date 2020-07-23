import metadataMap from '../../data/metadata_map.json';

const metadataFields = [
  'gender_id',
  'patient_status_id',
  'passage_id',
  'specimen_id',
  'sequencing_tech_id',
  'assembly_method_id',
  'comment_type_id',
];

const metadataFieldNiceNameMap = {
  gender_id: 'Gender',
  patient_status_id: 'Patient Status',
  passage_id: 'Passage',
  specimen_id: 'Specimen',
  sequencing_tech_id: 'Sequencing',
  assembly_method_id: 'Assembly',
  comment_type_id: 'Flag',
};

export function getMetadataFields() {
  return metadataFields;
}

export function getMetadataFieldNiceName(field) {
  return metadataFieldNiceNameMap[field];
}

export function countMetadataFields(caseData) {
  const metadataCounts = {};

  // Initialize fields
  metadataFields.forEach((field) => {
    metadataCounts[field] = {};
  });

  caseData.forEach((row) => {
    metadataFields.forEach((field) => {
      if (
        !Object.prototype.hasOwnProperty.call(metadataCounts[field], row[field])
      ) {
        metadataCounts[field][row[field]] = 0;
      }
      metadataCounts[field][row[field]] += 1;
    });
  });

  return metadataCounts;
}

export function getMetadataValueFromId(field, id) {
  return metadataMap[field.slice(0, -3)][id];
}
