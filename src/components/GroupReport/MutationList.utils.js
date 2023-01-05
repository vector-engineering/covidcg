import { DNA_OR_AA } from '../../constants/defs.json';
import { getAllGenes, getAllProteins } from '../../utils/gene_protein';
import { formatMutation } from '../../utils/mutationUtils';

const sortByPosThenAlt = function (a, b) {
  if (a.pos === b.pos) {
    return a.alt > b.alt;
  } else {
    return a.pos - b.pos;
  }
};

export const buildFeatureMatrix = ({
  reportActiveReference,
  reportGroupMutationFrequency,
  activeReportGroupType,
  reportGroupMutationType,
  selectedReportGroups,
  reportConsensusThreshold,
}) => {
  // Restructure so that we have it in a matrix-ish format
  /*
  {
    // gene or protein
    S: [
      {
        name: D614G
        ref,
        alt,
        pos: 614
        frequency: [0.1, 0.9, 0.5] // fractional frequencies per group
      },
      ...
    ],
    ...
  ]
  */

  // Abort if mutation frequency data empty
  if (Object.keys(reportGroupMutationFrequency).length === 0) {
    return {};
  }

  // Abort if mutation type data not fetched yet
  if (
    reportGroupMutationFrequency[activeReportGroupType] === undefined ||
    reportGroupMutationFrequency[activeReportGroupType][
      reportGroupMutationType
    ] === undefined
  ) {
    return {};
  }

  // Select group mutations from the selected groups
  reportGroupMutationFrequency = reportGroupMutationFrequency[
    activeReportGroupType
  ][reportGroupMutationType]['0'].filter((groupMutation) =>
    selectedReportGroups.includes(groupMutation.name)
  );

  const genes = getAllGenes(reportActiveReference);
  const proteins = getAllProteins(reportActiveReference);

  const features = reportGroupMutationType === 'protein_aa' ? proteins : genes;

  const result = {};

  features.forEach((feature) => {
    // Get all mutations for this gene, then sort by position/alt
    const groupFeatureMutations = reportGroupMutationFrequency
      .filter((groupMutation) => {
        if (reportGroupMutationType === 'dna') {
          // Include NT mutations in this gene if it is contained in
          // any of the gene's NT segments
          // (Most genes will have one segment)
          return feature.segments.some(
            (featureNTRange) =>
              groupMutation.pos >= featureNTRange[0] &&
              groupMutation.pos <= featureNTRange[1]
          );
        } else if (
          reportGroupMutationType === 'gene_aa' ||
          reportGroupMutationType == 'protein_aa'
        ) {
          return groupMutation.feature === feature.name;
        }
      })
      .sort(sortByPosThenAlt);

    // Make list of records to insert into master "matrix"
    const featureMutationRecords = groupFeatureMutations
      .slice()
      // Unique mutations
      .filter(
        (v, i, a) =>
          a.findIndex((element) => element.mutation_str === v.mutation_str) ===
          i
      )
      .map((featureMutation) => {
        // Find fractional frequencies for each group
        const freqs = [];
        selectedReportGroups.forEach((group) => {
          const matchingMutation = groupFeatureMutations.find(
            (mut) =>
              mut.mutation_str === featureMutation.mutation_str &&
              mut.name === group
          );
          // 0 if the mutation record isn't found
          freqs.push(
            matchingMutation === undefined ? 0 : matchingMutation.fraction
          );
        });

        return {
          // mutation_name: featureMutation.mutation_name,
          mutation_name: formatMutation(
            featureMutation.mutation_str,
            reportGroupMutationType === 'dna' ? DNA_OR_AA.DNA : DNA_OR_AA.AA
          ),
          pos: featureMutation.pos,
          ref: featureMutation.ref,
          alt: featureMutation.alt,
          frequency: freqs,
        };
      })
      // Filter out mutations that have all mutation frequencies below the threshold
      .filter((row) => {
        return row.frequency.some((freq) => freq > reportConsensusThreshold);
      });

    // console.log(featureMutationRecords);
    result[feature.name] = featureMutationRecords;
  });

  return result;
};

export const serializeFeatureMatrix = ({
  featureMatrix,
  selectedReportGroups,
}) => {
  let csvOut =
    'feature,mutation_name,ref,pos,alt,' +
    selectedReportGroups.join(',') +
    '\n';

  Object.keys(featureMatrix).forEach((featureName) => {
    const featureMutationRecords = featureMatrix[featureName];

    featureMutationRecords.forEach((mut) => {
      // Mutation metadata
      csvOut +=
        `${featureName},${mut.mutation_name},${mut.ref},${mut.pos},${mut.alt},` +
        mut.frequency.map((freq) => freq.toFixed(2)).join(',') +
        '\n';
    });
  });

  return csvOut;
};
