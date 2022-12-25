import {
  observable,
  action,
  toJS,
  //intercept, autorun
} from 'mobx';

import { coordsToText, residueCoordsToText } from '../utils/coordinates';
import { arrayEqual } from '../utils/func';
import { queryReferenceSequence } from '../utils/reference';
import { updateURLFromParams } from '../utils/updateQueryParam';

import {
  GEO_LEVELS,
  GROUP_MUTATION,
  DNA_OR_AA,
  COORDINATE_MODES,
  GROUPS,
} from '../constants/defs.json';
import { config } from '../config';

import { rootStoreInstance } from './rootStore';
import { configStore as initialConfigStore } from '../constants/initialValues';

export class ConfigStore {
  // Maintain a reference to the initial values
  // Initalize values to correct data types
  initialValues = {};

  @observable groupKey = '';
  @observable dnaOrAa = '';

  @observable selectedReference = '';

  @observable selectedGene = {};
  @observable selectedProtein = {};
  @observable selectedPrimers = [];

  @observable customCoordinates = [[]];
  @observable customSequences = [];
  @observable residueCoordinates = [[]];
  @observable coordinateMode = '';

  @observable startDate = new Date();
  @observable endDate = new Date();

  @observable submStartDate = '';
  @observable submEndDate = '';

  @observable selectedGroupFields = {};
  @observable selectedLocationNodes = [];

  @observable hoverGroup = null;
  @observable selectedGroups = [];

  @observable selectedMetadataFields = {};
  @observable ageRange = [];

  @observable hoverLocation = null;
  @observable focusedLocations = [];

  constructor() {}

  init() {
    this.initialValues = initialConfigStore;

    Object.keys(this.initialValues).forEach((key) => {
      this[key] = this.initialValues[key];
    });
  }

  @action
  resetValues = (values) => {
    Object.keys(this.initialValues).forEach((key) => {
      if (key in values) {
        this[key] = values[key];
      } else {
        this[key] = this.initialValues[key];
      }

      // Special actions for some keys
      if (key === 'selectedLocationNodes') {
        rootStoreInstance.locationDataStore.setSelectedNodes(values[key]);
      }
    });

    // Manually set residue coordinates if they weren't specified
    if (
      'selectedGene' in values &&
      !('residueCoordinates' in values) &&
      this.selectedGene.name !== 'All Genes'
    ) {
      this.residueCoordinates = [[1, this.selectedGene.len_aa]];
    } else if (
      'selectedProtein' in values &&
      !('residueCoordinates' in values) &&
      this.selectedProtein.name !== 'All Proteins'
    ) {
      this.residueCoordinates = [[1, this.selectedProtein.len_aa]];
    }

    // Trigger data re-run
    rootStoreInstance.dataStore.fetchData();
  };

  @action
  applyPendingChanges = (pending, fetch = true) => {
    // Clear selected groups/locations
    this.hoverGroup = this.initialValues.hoverGroup;
    this.selectedGroups = this.initialValues.selectedGroups;
    this.hoverLocation = this.initialValues.hoverLocation;
    this.focusedLocations = this.initialValues.focusedLocations;

    // Overwrite any of our fields here with the pending ones
    Object.keys(pending).forEach((field) => {
      this[field] = pending[field];
    });

    // Update the location node tree with our new selection
    rootStoreInstance.locationDataStore.setSelectedNodes(
      this.selectedLocationNodes
    );

    // Get the new data from the server
    if (fetch) {
      rootStoreInstance.dataStore.fetchData();
    }

    this.updateURL(pending);
  };

  /*
   * Serialize store state into URL params
   */
  updateURL = (pending) => {
    const urlParams = rootStoreInstance.urlMonitor.urlParams;

    Object.keys(pending).forEach((field) => {
      if (field === 'selectedGene' || field === 'selectedProtein') {
        // Handle fields that return objects
        urlParams.set(field, pending[field].name);
      } else if (
        field === 'selectedMetadataFields' ||
        field === 'selectedLocationNodes' ||
        field === 'ageRange'
      ) {
        // Ignore Metadata fields
        // selectedLocationNodes is displayed as node.level=node.value in the URL
        // ageRange is not currently being used
        return;
      } else if (field.includes('valid')) {
        // Ignore boolean flags
        return;
      } else if (field === 'selectedPrimers') {
        pending[field].forEach((primer) => {
          if (urlParams.has(field)) {
            urlParams.append(field, primer.Institution + '_' + primer.Name);
          } else {
            urlParams.set(field, primer.Institution + '_' + primer.Name);
          }
        });
      } else if (field === 'selectedGroupFields') {
        // If selectedGroupFields just got emptied - make sure to
        // remove any lingering group keys in the URL
        Object.keys(config.group_cols).forEach((groupKey) => {
          if (!Object.keys(pending[field]).includes(groupKey)) {
            urlParams.delete(groupKey);
          }
        });

        Object.keys(pending[field]).forEach((groupKey) => {
          urlParams.delete(groupKey);
          pending[field][groupKey].forEach((group) => {
            urlParams.append(groupKey, group);
          });
        });
      } else if (field === 'customCoordinates') {
        urlParams.set(field, coordsToText(pending[field]));
      } else if (field === 'residueCoordinates') {
        urlParams.set(field, residueCoordsToText(pending[field]));
      } else {
        urlParams.set(field, String(pending[field]));
      }

      if (
        JSON.stringify(pending[field]) ===
        JSON.stringify(this.initialValues[field])
      ) {
        // Only display non-default fields in the url
        urlParams.delete(field);
      }
    });

    // Show only relevant coordinate info
    const coordinateMode = urlParams.get('coordinateMode');
    switch (coordinateMode) {
      case 'protein':
        urlParams.delete('selectedGene');
        urlParams.delete('selectedPrimers');
        urlParams.delete('customCoordinates');
        urlParams.delete('customSequences');
        break;
      // Gene is the default coordinateMode and may have been deleted so check for null
      case 'gene':
      case null:
      case undefined:
        urlParams.delete('selectedProtein');
        urlParams.delete('selectedPrimers');
        urlParams.delete('customCoordinates');
        urlParams.delete('customSequences');
        break;
      case 'primer':
        urlParams.delete('selectedProtein');
        urlParams.delete('selectedGene');
        urlParams.delete('customCoordinates');
        urlParams.delete('customSequences');
        break;
      case 'custom':
        urlParams.delete('selectedProtein');
        urlParams.delete('selectedPrimers');
        urlParams.delete('selectedGene');
        urlParams.delete('customSequences');
        break;
      case 'sequence':
        urlParams.delete('selectedProtein');
        urlParams.delete('selectedPrimers');
        urlParams.delete('selectedGene');
        urlParams.delete('customCoordinates');
        break;
    }

    if (this.groupKey !== GROUP_MUTATION) {
      urlParams.delete('residueCoordinates');
      urlParams.delete('selectedGene');
      urlParams.delete('selectedProtein');
      urlParams.delete('selectedPrimers');
      urlParams.delete('customCoordinates');
      urlParams.delete('customSequences');
    }

    // Update the location URL params
    Object.values(GEO_LEVELS).forEach((level) => {
      urlParams.delete(level);
    });

    if ('selectedLocationNodes' in pending) {
      pending.selectedLocationNodes.forEach((node) => {
        if (urlParams.has(String(node.level))) {
          urlParams.append(String(node.level), String(node.value_txt));
        } else {
          urlParams.set(String(node.level), String(node.value_txt));
        }
      });
    }

    // Update URL
    updateURLFromParams(urlParams);

    rootStoreInstance.urlMonitor.urlParams = urlParams;
  };

  getMutationType() {
    if (this.dnaOrAa === DNA_OR_AA.DNA) {
      return 'dna';
    } else if (this.dnaOrAa === DNA_OR_AA.AA) {
      if (this.coordinateMode === COORDINATE_MODES.COORD_GENE) {
        return 'gene_aa';
      } else if (this.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
        return 'protein_aa';
      }
    }
    return undefined;
  }

  // Get a pretty name for the group
  getGroupLabel(groupKey = null, dnaOrAa = null) {
    // Default to using store attribute, if no explicit groupKey
    // is provided
    if (groupKey === null) {
      groupKey = this.groupKey;
    }
    if (dnaOrAa === null) {
      dnaOrAa = this.dnaOrAa;
    }

    if (Object.keys(config.group_cols).includes(groupKey)) {
      return config.group_cols[groupKey].title;
    } else if (groupKey === GROUP_MUTATION) {
      if (dnaOrAa === DNA_OR_AA.DNA) {
        return 'NT Mutation';
      } else {
        return 'AA Mutation';
      }
    }
  }

  getSelectedLocations() {
    const res = {
      region: [],
      country: [],
      division: [],
      location: [],
    };
    this.selectedLocationNodes.forEach((node) => {
      res[node.level].push(node.value);
    });
    return res;
  }

  getCoordinateRanges() {
    // Set the coordinate range based off the coordinate mode
    if (
      this.coordinateMode === COORDINATE_MODES.COORD_GENE ||
      this.coordinateMode === COORDINATE_MODES.COORD_PROTEIN
    ) {
      const feature =
        this.coordinateMode === COORDINATE_MODES.COORD_GENE
          ? this.selectedGene
          : this.selectedProtein;

      // Return ranges if All Genes
      // if (feature.name === 'All Genes') {
      //   return feature.ranges;
      // }
      // Disable residue indices for non-protein-coding genes
      if (!feature.protein_coding) {
        return feature.segments.slice().map((segment) => {
          return [feature.segment, segment[0], segment[1]];
        });
      }
      const coordinateRanges = [];
      this.residueCoordinates.forEach((range) => {
        // Make a deep copy of the current range
        const curRange = range.slice();

        if (this.dnaOrAa === DNA_OR_AA.DNA) {
          for (let i = 0; i < feature.aa_ranges.length; i++) {
            const curAARange = feature.aa_ranges[i];
            const curNTRange = feature.segments[i];
            if (
              (curRange[0] >= curAARange[0] && curRange[0] <= curAARange[1]) ||
              (curRange[0] <= curAARange[0] && curRange[1] >= curAARange[0])
            ) {
              coordinateRanges.push([
                feature.segment,
                curNTRange[0] + (curRange[0] - curAARange[0]) * 3,
                curNTRange[0] -
                  1 +
                  Math.min(curRange[1] - curAARange[0] + 1, curAARange[1]) * 3,
              ]);
              // Push the beginning of the current range to the end of
              // the current AA range of the gene
              if (curAARange[1] < curRange[1]) {
                curRange[0] = curAARange[1] + 1;
              }
            }
          }
        } else {
          coordinateRanges.push([feature.segment, curRange[0], curRange[1]]);
        }
      });
      return coordinateRanges;
    } else if (this.coordinateMode === COORDINATE_MODES.COORD_PRIMER) {
      return this.selectedPrimers.map((primer) => {
        return [primer.Start, primer.End];
      });
    } else if (this.coordinateMode === COORDINATE_MODES.COORD_CUSTOM) {
      return toJS(this.customCoordinates);
    } else if (this.coordinateMode === COORDINATE_MODES.COORD_SEQUENCE) {
      return this.customSequences.map((seq) => {
        return queryReferenceSequence(seq, this.selectedReference);
      });
    }
  }

  getSelectedMetadataFields() {
    const selectedMetadataFields = toJS(this.selectedMetadataFields);
    Object.keys(selectedMetadataFields).forEach((metadataField) => {
      selectedMetadataFields[metadataField] = selectedMetadataFields[
        metadataField
      ].map((item) => {
        return parseInt(item.value);
      });
    });
    return selectedMetadataFields;
  }

  @action
  updateHoverGroup = (group) => {
    if (group === this.hoverGroup) {
      return;
    } else if (group === GROUPS.NONE_GROUP || GROUPS.ALL_OTHER_GROUP) {
      // Ignore for some special groups
      return;
    } else if (group === null) {
      this.hoverGroup = null;
    } else {
      this.hoverGroup = group;
    }
  };

  @action
  updateSelectedGroups = (groups) => {
    // First check to see that it's different. If not,
    // skip the update
    // This matters because JS arrays are passed by reference,
    // and any listeners which see a change to this reference
    // will update themselves when the reference changes,
    // even if the underlying data does not
    // groups are in the structure [{ group: group }]
    if (
      arrayEqual(
        groups.map((group) => group.group),
        this.selectedGroups.map((group) => group.group)
      )
    ) {
      return;
    }

    this.selectedGroups = groups;
    if (this.groupKey === GROUP_MUTATION) {
      rootStoreInstance.dataStore.processSelectedMutations();
    }
  };

  getSelectedGroupIds() {
    const { dnaMutationMap, geneAaMutationMap, proteinAaMutationMap } =
      rootStoreInstance.mutationDataStore;

    let selectedGroupIds;
    if (this.dnaOrAa === DNA_OR_AA.DNA) {
      selectedGroupIds = this.selectedGroups
        .map((item) => dnaMutationMap[item.group])
        .map((mutationId) =>
          mutationId === undefined ? -1 : parseInt(mutationId)
        );
    } else if (this.dnaOrAa === DNA_OR_AA.AA) {
      if (this.coordinateMode === COORDINATE_MODES.COORD_GENE) {
        selectedGroupIds = this.selectedGroups
          .map((item) => geneAaMutationMap[item.group])
          .map((mutationId) =>
            mutationId === undefined ? -1 : parseInt(mutationId)
          );
      } else if (this.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
        selectedGroupIds = this.selectedGroups
          .map((item) => proteinAaMutationMap[item.group])
          .map((mutationId) =>
            mutationId === undefined ? -1 : parseInt(mutationId)
          );
      }
    }
    // Array to Set
    selectedGroupIds = new Set(selectedGroupIds);

    return selectedGroupIds;
  }

  getIntToMutationMap() {
    const {
      intToDnaMutationMap,
      intToGeneAaMutationMap,
      intToProteinAaMutationMap,
    } = rootStoreInstance.mutationDataStore;

    if (this.dnaOrAa === DNA_OR_AA.DNA) {
      return intToDnaMutationMap;
    } else if (this.dnaOrAa === DNA_OR_AA.AA) {
      if (this.coordinateMode === COORDINATE_MODES.COORD_GENE) {
        return intToGeneAaMutationMap;
      } else if (this.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
        return intToProteinAaMutationMap;
      }
    }
  }

  getMutationToIntMap() {
    const { dnaMutationMap, geneAaMutationMap, proteinAaMutationMap } =
      rootStoreInstance.mutationDataStore;

    if (this.dnaOrAa === DNA_OR_AA.DNA) {
      return dnaMutationMap;
    } else if (this.dnaOrAa === DNA_OR_AA.AA) {
      if (this.coordinateMode === COORDINATE_MODES.COORD_GENE) {
        return geneAaMutationMap;
      } else if (this.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
        return proteinAaMutationMap;
      }
    }
  }

  @action
  updateHoverLocation = (location) => {
    this.hoverLocation = location;
  };

  @action
  updateFocusedLocations = (locations) => {
    if (
      arrayEqual(
        locations.map((location) => location.location),
        this.focusedLocations.map((location) => location.location)
      )
    ) {
      return;
    }
    this.focusedLocations = locations;
  };
}
