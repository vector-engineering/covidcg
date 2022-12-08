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
  activeReference,
  groupMutationFrequency,
  activeGroupType,
  groupMutationType,
  selectedGroups,
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
  if (Object.keys(groupMutationFrequency).length === 0) {
    return {};
  }

  // Select group mutations from the selected groups
  groupMutationFrequency = groupMutationFrequency[activeGroupType][
    groupMutationType
  ]['0'].filter((groupMutation) => selectedGroups.includes(groupMutation.name));

  const genes = getAllGenes(activeReference);
  const proteins = getAllProteins(activeReference);

  const features = groupMutationType === 'protein_aa' ? proteins : genes;

  const result = {};

  features.forEach((feature) => {
    // Get all mutations for this gene, then sort by position/alt
    const groupFeatureMutations = groupMutationFrequency
      .filter((groupMutation) => {
        if (groupMutationType === 'dna') {
          // Include NT mutations in this gene if it is contained in
          // any of the gene's NT segments
          // (Most genes will have one segment)
          return feature.segments.some(
            (featureNTRange) =>
              groupMutation.pos >= featureNTRange[0] &&
              groupMutation.pos <= featureNTRange[1]
          );
        } else if (
          groupMutationType === 'gene_aa' ||
          groupMutationType == 'protein_aa'
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
        selectedGroups.forEach((group) => {
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
            groupMutationType === 'dna' ? DNA_OR_AA.DNA : DNA_OR_AA.AA
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

export const serializeFeatureMatrix = ({ featureMatrix, selectedGroups }) => {
  let csvOut =
    'feature,mutation_name,ref,pos,alt,' + selectedGroups.join(',') + '\n';

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
