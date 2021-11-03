import { getGene, getProtein } from '../utils/gene_protein';
import { intToISO, ISOToInt } from '../utils/date';

import { DNA_OR_AA, COORDINATE_MODES } from '../constants/defs.json';

const today = intToISO(new Date().getTime());
const lastNDays = 30; // By default, show only the last 1 month

export default function values() {
  return {
    configStore: {
      groupKey: 'snv',
      dnaOrAa: DNA_OR_AA.AA,

      // Select the Spike gene and nsp13 protein by default
      selectedGene: getGene('S'),
      selectedProtein: getProtein('nsp12 - RdRp'),
      selectedPrimers: [],
      customCoordinates: [[8000, 12000]],
      customSequences: ['GACCCCAAAATCAGCGAAAT'],
      residueCoordinates: [[1, getGene('S').len_aa]],

      // Selecting the gene as the coordinate range by default
      coordinateMode: COORDINATE_MODES.COORD_GENE,

      // days * (24 hours/day) * (60 min/hour) * (60 s/min) * (1000 ms/s)
      startDate: intToISO(ISOToInt(today) - lastNDays * 24 * 60 * 60 * 1000),
      endDate: today,

      submStartDate: '',
      submEndDate: '',

      selectedLocationNodes: [],

      hoverGroup: null,
      selectedGroups: [],

      // Metadata filtering
      selectedMetadataFields: {},
      ageRange: [null, null],

      // Location tab
      hoverLocation: null,
      focusedLocations: [],
    },
  };
}
