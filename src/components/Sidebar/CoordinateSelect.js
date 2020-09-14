import React, { useState, useEffect, useMemo } from 'react';
// import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import _ from 'underscore';

import ReactTooltip from 'react-tooltip';
import ExternalLink from '../Common/ExternalLink';
import QuestionButton from '../Buttons/QuestionButton';
import {
  PrimerSelect,
  ModeLabel,
  SelectContainer,
  ModeSelectForm,
  ModeRadioVertical,
  SelectForm,
  UpdatePrimersButton,
  UpdateButton,
  CoordForm,
  PrimerSelectContainer,
  ValidationInput,
  InvalidText,
  RangesText,
  DomainSelectForm,
} from './CoordinateSelect.styles';

import { getAllGenes, getAllProteins } from '../../utils/gene_protein';
import {
  getPrimerSelectTree,
  getPrimerByName,
  getPrimersByGroup,
} from '../../utils/primer';
import { referenceSequenceIncludes } from '../../utils/reference';

import { DNA_OR_AA, COORDINATE_MODES } from '../../constants/config';

const genes = getAllGenes();
const proteins = getAllProteins();

const CoordinateSelect = observer(() => {
  const { configStore } = useStores();

  // Create option elements

  // GENE
  let geneOptionElements = [];
  genes.forEach((gene) => {
    geneOptionElements.push(
      <option
        key={gene.gene}
        value={gene.gene}
        disabled={
          gene.protein_coding === 0 && configStore.dnaOrAa === DNA_OR_AA.AA
        }
      >
        {gene.gene}&nbsp;&nbsp;({gene.segments})
      </option>
    );
  });

  // GENE DOMAINS
  let geneDomainOptionElements = {};
  genes.forEach((gene) => {
    geneDomainOptionElements[gene.gene] = [
      <option
        key={`${gene.gene}-default`}
        value={`${gene.gene}-default`}
        disabled={true}
      >
        - select an option -
      </option>,
      <option key={`${gene.gene}-all`} value={`${gene.gene}-all`}>
        Entire {gene.gene} Gene (1..{gene.len_aa})
      </option>,
    ];
    gene.domains.forEach((domain) => {
      geneDomainOptionElements[gene.gene].push(
        <option key={`${gene.gene}-${domain.name}`} value={domain.name}>
          {domain.name}&nbsp;&nbsp;(
          {domain.ranges.map((range) => range.join('..')).join(';')})
        </option>
      );
    });
  });

  // PROTEIN
  let proteinOptionElements = [];
  proteins.forEach((protein) => {
    proteinOptionElements.push(
      <option key={protein.protein} value={protein.protein}>
        {protein.protein}&nbsp;&nbsp;({protein.segments})
      </option>
    );
  });

  // PROTEIN DOMAINS
  let proteinDomainOptionElements = {};
  proteins.forEach((protein) => {
    proteinDomainOptionElements[protein.protein] = [
      <option
        key={`${protein.protein}-default`}
        value={`${protein.protein}-default`}
        disabled={true}
      >
        {' '}
        - select an option -
      </option>,
      <option key={`${protein.protein}-all`} value={`${protein.protein}-all`}>
        Entire {protein.protein} Protein (1..{protein.len_aa})
      </option>,
    ];
    protein.domains.forEach((domain) => {
      proteinDomainOptionElements[protein.protein].push(
        <option key={`${protein.protein}-${domain.name}`} value={domain.name}>
          {domain.name}&nbsp;&nbsp;(
          {domain.ranges.map((range) => range.join('..')).join(';')})
        </option>
      );
    });
  });

  const [state, setState] = useState({
    primerTreeData: Object.assign(getPrimerSelectTree()),
    selectedPrimers: [],
    primersChanged: false,

    customCoordText: configStore.customCoordinates
      .map((range) => range.join('..'))
      .join(';'),
    validCustomCoords: true,
    customCoordinatesChanged: false,

    customSequences: configStore.customSequences.join(';'),
    validCustomSequences: true,
    customSequencesChanged: false,

    residueCoordsText: configStore.residueCoordinates
      .map((range) => range.join('..'))
      .join(';'),
    validResidueCoords: true,
    residueCoordsChanged: false,
  });

  // Disable "All Genes" and "All Proteins" option
  // when in AA mode and non-SNP grouping
  // useEffect(() => {
  //   let _geneOptionElements = state.geneOptionElements;
  //   let _proteinOptionElements = state.proteinOptionElements;
  //   if (configStore.groupKey !== GROUP_KEYS.GROUP_SNV && configStore.dnaOrAa === DNA_OR_AA.AA) {
  //   }
  // }, [configStore.groupKey, configStore.dnaOrAa]);

  // Update custom coordinates from the store
  useEffect(() => {
    setState({
      ...state,
      customCoordText: configStore.customCoordinates
        .map((range) => range.join('..'))
        .join(';'),
      customCoordinatesChanged: false,
      validCustomCoords: true,
    });
  }, [configStore.customCoordinates]);

  // Update custom sequences from the store
  useEffect(() => {
    setState({
      ...state,
      customSequences: configStore.customSequences.join(';'),
      customSequencesChanged: false,
      validCustomSequences: true,
    });
  }, [configStore.customSequences]);

  const changeCoordinateMode = ({
    coordinateMode,
    selectedGene,
    selectedProtein,
    residueCoordinates,
    selectedPrimers,
    customCoordinates,
    customSequences,
  }) => {
    configStore.changeCoordinateMode({
      coordinateMode:
        coordinateMode === undefined
          ? configStore.coordinateMode
          : coordinateMode,
      selectedGene:
        selectedGene === undefined
          ? configStore.selectedGene.gene
          : selectedGene,
      selectedProtein:
        selectedProtein === undefined
          ? configStore.selectedProtein.protein
          : selectedProtein,
      residueCoordinates:
        residueCoordinates === undefined
          ? configStore.residueCoordinates
          : residueCoordinates,
      selectedPrimers:
        selectedPrimers === undefined
          ? configStore.selectedPrimers
          : selectedPrimers,
      customCoordinates:
        customCoordinates === undefined
          ? configStore.customCoordinates
          : customCoordinates,
      customSequences:
        customSequences === undefined
          ? configStore.customSequences
          : customSequences,
    });
  };

  const handleModeChange = (event) => {
    changeCoordinateMode({ coordinateMode: event.target.value });
  };

  const handleGeneChange = (event) => {
    changeCoordinateMode({
      coordinateMode: COORDINATE_MODES.COORD_GENE,
      selectedGene: event.target.value,
    });
  };

  const handleProteinChange = (event) => {
    changeCoordinateMode({
      coordinateMode: COORDINATE_MODES.COORD_PROTEIN,
      selectedProtein: event.target.value,
    });
  };

  // Use a regex to match numbers, since just because JS
  // can parse an integer, doesn't mean it should...
  const numPattern = /^([0-9]+)$/;

  const handleResidueCoordsChange = (event) => {
    // Serialize custom coordinates of the store
    const storeResidueCoords = configStore.residueCoordinates
      .map((range) => range.join('..'))
      .join(';');

    // Parse current custom coordinates
    const curResidueCoords = event.target.value
      .split(';')
      .map((range) => range.split('..'));

    // Check that these are valid
    const validResidueCoords = !curResidueCoords.some((range) => {
      // Return true if invalid
      return (
        range.length !== 2 ||
        numPattern.exec(range[0]) === null ||
        numPattern.exec(range[1]) === null ||
        parseInt(range[0]) > parseInt(range[1])
      );
    });

    setState({
      ...state,
      residueCoordsText: event.target.value,
      validResidueCoords,
      residueCoordsChanged: storeResidueCoords !== event.target.value,
    });
  };

  const handleResidueCoordsSubmit = (event) => {
    event.preventDefault();
    changeCoordinateMode({
      residueCoordinates: state.residueCoordsText
        .split(';')
        .map((range) => range.split('..').map((coord) => parseInt(coord))),
    });
  };

  // Use the selected domain to fill in the residue coordinates input
  const handleGeneDomainChange = (event) => {
    const domainName = event.target.value;
    let newResidueCoordsText;

    if (event.target.value === configStore.selectedGene.gene + '-all') {
      newResidueCoordsText = `1..${configStore.selectedGene.len_aa}`;
    } else {
      const domainObj = _.findWhere(configStore.selectedGene.domains, {
        name: domainName,
      });

      newResidueCoordsText = domainObj.ranges
        .map((range) => range.join('..'))
        .join(';');
    }

    const residueCoordsChanged =
      state.residueCoordsText !== newResidueCoordsText;

    setState({
      ...state,
      residueCoordsText: newResidueCoordsText,
      validResidueCoords: true,
      residueCoordsChanged,
    });
  };

  const handleProteinDomainChange = (event) => {
    const domainName = event.target.value;
    let newResidueCoordsText;

    if (event.target.value === configStore.selectedProtein.protein + '-all') {
      newResidueCoordsText = `1..${configStore.selectedProtein.len_aa}`;
    } else {
      const domainObj = _.findWhere(configStore.selectedProtein.domains, {
        name: domainName,
      });

      newResidueCoordsText = domainObj.ranges
        .map((range) => range.join('..'))
        .join(';');
    }
    const residueCoordsChanged =
      state.residueCoordsText !== newResidueCoordsText;

    setState({
      ...state,
      residueCoordsText: newResidueCoordsText,
      validResidueCoords: true,
      residueCoordsChanged,
    });
  };

  // Update residue coordinates from store
  useEffect(() => {
    setState({
      ...state,
      residueCoordsText: configStore.residueCoordinates
        .map((range) => range.join('..'))
        .join(';'),
      validResidueCoords: true,
      residueCoordsChanged: false,
    });
  }, [configStore.residueCoordinates]);

  const handleCustomCoordChange = (event) => {
    // Serialize custom coordinates of the store
    const storeCustomCoords = configStore.customCoordinates
      .map((range) => range.join('..'))
      .join(';');

    // Parse current custom coordinates
    const curCustomCoords = event.target.value
      .split(';')
      .map((range) => range.split('..'));
    // Check that these are valid
    const validCustomCoords = !curCustomCoords.some((range) => {
      // Return true if invalid
      return (
        range.length !== 2 ||
        numPattern.exec(range[0]) === null ||
        numPattern.exec(range[1]) === null ||
        parseInt(range[0]) > parseInt(range[1])
      );
    });

    setState({
      ...state,
      customCoordText: event.target.value,
      validCustomCoords,
      customCoordinatesChanged: storeCustomCoords !== event.target.value,
    });
  };

  const handleCustomCoordSubmit = (event) => {
    event.preventDefault();
    // Change to custom mode implicitly
    changeCoordinateMode({
      coordinateMode: COORDINATE_MODES.COORD_CUSTOM,
      customCoordinates: state.customCoordText
        .split(';')
        .map((range) => range.split('..').map((coord) => parseInt(coord))),
    });
  };

  const handleCustomSequencesChange = (event) => {
    const curText = event.target.value.toUpperCase();
    const sequences = curText.split(';');
    // Check that all bases are valid and that
    // the reference sequence includes the sequence
    // TODO: support for denegerate bases
    const validBases = ['A', 'T', 'C', 'G'];
    const validCustomSequences = !sequences.some((seq) => {
      // Fails if any conditions are met
      return (
        seq.length === 0 ||
        Array.from(seq).some((base) => !validBases.includes(base)) ||
        !referenceSequenceIncludes(seq)
      );
    });
    const customSequencesChanged =
      curText !== configStore.customSequences.join(';');

    setState({
      ...state,
      customSequences: curText,
      validCustomSequences,
      customSequencesChanged,
    });
  };

  const handleCustomSequencesSubmit = (event) => {
    event.preventDefault();
    // Change to sequences mode implicitly
    changeCoordinateMode({
      coordinateMode: COORDINATE_MODES.COORD_SEQUENCE,
      customSequences: state.customSequences.split(';'),
    });
  };

  const checkPrimersChanged = (selectedPrimers) => {
    // Is the current configStore.selectedPrimers the same as the current selection?
    let changed = selectedPrimers.length !== configStore.selectedPrimers.length;
    if (!changed && selectedPrimers.length > 0) {
      // Run through once - both selectedPrimers and the configStore version are sorted
      for (let i = 0; i < selectedPrimers.length; i++) {
        if (!_.isEqual(selectedPrimers[i], configStore.selectedPrimers[i])) {
          changed = true;
          break;
        }
      }
    }
    return changed;
  };

  useEffect(() => {
    // Check all selected primers
    // Make a deep copy of the primer tree data - so we trigger an update
    // in the memoized primer tree element
    const primerTreeData = state.primerTreeData.slice();

    // Recursively go through and deselect everything
    const traverseAndDeselect = (node) => {
      node.checked = false;
      if ('children' in node) {
        node.children.forEach((child) => {
          traverseAndDeselect(child);
        });
      }
    };
    primerTreeData.forEach((node) => {
      traverseAndDeselect(node);
    });

    const selectedPrimers = configStore.selectedPrimers;
    selectedPrimers.forEach((primer) => {
      const institutionNode = _.findWhere(primerTreeData, {
        value: primer.Institution,
      });
      const primerNode = _.findWhere(institutionNode.children, {
        value: primer.Name,
      });
      primerNode.checked = true;
    });

    setState({
      ...state,
      primerTreeData,
      selectedPrimers,
      primersChanged: false,
    });
  }, [configStore.selectedPrimers]);

  const onPrimerSelect = (currentNode, selectedNodes) => {
    // console.log(currentNode);
    // console.log(selectedNodes);

    let selectedPrimers = [];
    selectedNodes.forEach((node) => {
      if (node.level === 'group') {
        selectedPrimers = selectedPrimers.concat(getPrimersByGroup(node.value));
      } else if (node.level === 'individual') {
        selectedPrimers = selectedPrimers.concat(getPrimerByName(node.value));
      }
    });

    // Sort by Institution then Name
    selectedPrimers = _.sortBy(selectedPrimers, (primer) => {
      return primer.Institution.concat('-', primer.Name);
    });

    setState({
      ...state,
      selectedPrimers: selectedPrimers,
      primersChanged: checkPrimersChanged(selectedPrimers),
    });
  };

  const updatePrimerSelection = (event) => {
    event.preventDefault();
    // In addition to updating the selection, also
    // switch to primer mode here implicitly
    changeCoordinateMode({
      coordinateMode: COORDINATE_MODES.COORD_PRIMER,
      selectedPrimers: state.selectedPrimers,
    });
  };

  // This component needs to be in a memoized function
  // since it manages its own local state. It should never be re-rendered
  // forcefully
  const primerDropdown = useMemo(() => {
    return (
      <PrimerSelect
        data={state.primerTreeData}
        className="primer-dropdown-tree-select"
        clearSearchOnChange={false}
        keepTreeOnSearch={true}
        keepChildrenOnSearch={true}
        showPartiallySelected={true}
        inlineSearchInput={true}
        texts={{
          placeholder: 'Search...',
          noMatches: 'No matches found',
        }}
        onChange={onPrimerSelect}
      />
    );
  }, [state.primerTreeData]);

  return (
    <SelectContainer>
      <ModeSelectForm>
        {/* GENE SELECT */}
        <ModeRadioVertical>
          <ModeLabel>
            <input
              className="radio-input"
              type="radio"
              value={COORDINATE_MODES.COORD_GENE}
              checked={
                configStore.coordinateMode === COORDINATE_MODES.COORD_GENE
              }
              onChange={handleModeChange}
            />
            <span className="option-text">Gene</span>
            <SelectForm>
              <select
                value={configStore.selectedGene.gene}
                onChange={handleGeneChange}
              >
                <option
                  key="All Genes"
                  value="All Genes"
                  disabled={configStore.dnaOrAa === DNA_OR_AA.AA}
                >
                  All Genes
                </option>
                {geneOptionElements}
              </select>
            </SelectForm>
          </ModeLabel>
          {configStore.coordinateMode === COORDINATE_MODES.COORD_GENE &&
            configStore.selectedGene.gene !== 'All Genes' && (
              <>
                <CoordForm>
                  <span className="coord-prefix">Residue indices:</span>
                  <input
                    type="text"
                    value={state.residueCoordsText}
                    onChange={handleResidueCoordsChange}
                  />
                  <ReactTooltip
                    className="filter-sidebar-tooltip"
                    id="gene-residue-index-tooltip"
                    type="light"
                    effect="solid"
                    border={true}
                    borderColor="#888"
                  />
                  <QuestionButton
                    data-tip='<p>Coordinates are in the form "start..end". Multiple ranges can be separated with ";"</p><p>i.e., "100..300;500..550"</p><p>Coordinates are relative to the gene ORF</p>'
                    data-html="true"
                    data-for="gene-residue-index-tooltip"
                  />
                </CoordForm>
                <DomainSelectForm>
                  <span>Domain:</span>
                  <select
                    value={`${configStore.selectedGene.gene}-default`}
                    onChange={handleGeneDomainChange}
                  >
                    {geneDomainOptionElements[configStore.selectedGene.gene]}
                  </select>
                  <QuestionButton
                    data-tip='<p>Coordinates relative to the gene ORF, and are in the form "start..end".</p><p>Selecting a domain will replace the range(s) to the residue indices input</p>'
                    data-html="true"
                    data-for="gene-residue-index-tooltip"
                  />
                </DomainSelectForm>

                <UpdateButton
                  show={state.residueCoordsChanged}
                  disabled={!state.validResidueCoords}
                  onClick={handleResidueCoordsSubmit}
                  style={{ marginTop: 5 }}
                >
                  Confirm
                </UpdateButton>
              </>
            )}
        </ModeRadioVertical>

        {/* PROTEIN SELECT */}
        <ModeRadioVertical>
          <ModeLabel>
            <input
              className="radio-input"
              type="radio"
              value={COORDINATE_MODES.COORD_PROTEIN}
              checked={
                configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN
              }
              onChange={handleModeChange}
            />
            <span className="option-text">Protein</span>
            <SelectForm>
              <select
                value={configStore.selectedProtein.protein}
                onChange={handleProteinChange}
              >
                <option
                  key="All Proteins"
                  value="All Proteins"
                  disabled={configStore.dnaOrAa === DNA_OR_AA.AA}
                >
                  All Proteins
                </option>
                {proteinOptionElements}
              </select>
            </SelectForm>
          </ModeLabel>
          {configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN &&
            configStore.selectedProtein.protein !== 'All Proteins' && (
              <>
                <CoordForm>
                  <span className="coord-prefix">Residue indices:</span>
                  <input
                    type="text"
                    value={state.residueCoordsText}
                    onChange={handleResidueCoordsChange}
                  />
                  <ReactTooltip
                    className="filter-sidebar-tooltip"
                    id="protein-residue-index-tooltip"
                    type="light"
                    effect="solid"
                    border={true}
                    borderColor="#888"
                  />
                  <QuestionButton
                    data-tip='<p>Coordinates are in the form "start..end". Multiple ranges can be separated with ";"</p><p>i.e., "100..300;500..550"</p><p>Coordinates are relative to the protein ORF</p>'
                    data-html="true"
                    data-for="protein-residue-index-tooltip"
                  />
                </CoordForm>
                <DomainSelectForm>
                  <span>Domain:</span>
                  <select
                    value={`${configStore.selectedProtein.protein}-default`}
                    onChange={handleProteinDomainChange}
                  >
                    {
                      proteinDomainOptionElements[
                        configStore.selectedProtein.protein
                      ]
                    }
                  </select>
                  <QuestionButton
                    data-tip='<p>Coordinates relative to the protein ORF, and are in the form "start..end".</p><p>Selecting a domain will replace the range(s) to the residue indices input</p>'
                    data-html="true"
                    data-for="gene-residue-index-tooltip"
                  />
                </DomainSelectForm>
                <UpdateButton
                  show={state.residueCoordsChanged}
                  disabled={!state.validResidueCoords}
                  onClick={handleResidueCoordsSubmit}
                  style={{ marginTop: 5 }}
                >
                  Confirm
                </UpdateButton>
              </>
            )}
        </ModeRadioVertical>

        {/* PRIMER/PROBE SELECT */}
        <ModeRadioVertical>
          <ModeLabel>
            <input
              className="radio-input"
              type="radio"
              value={COORDINATE_MODES.COORD_PRIMER}
              checked={
                configStore.coordinateMode === COORDINATE_MODES.COORD_PRIMER
              }
              onChange={handleModeChange}
            />
            <span className="select-text">Primers/Probes</span>
            {configStore.coordinateMode !== COORDINATE_MODES.COORD_PRIMER && (
              <span className="hint-text">Select to show options</span>
            )}
            {configStore.coordinateMode === COORDINATE_MODES.COORD_PRIMER && (
              <UpdatePrimersButton
                show={state.primersChanged}
                onClick={updatePrimerSelection}
              >
                Update Selection
              </UpdatePrimersButton>
            )}
          </ModeLabel>
          {configStore.coordinateMode === COORDINATE_MODES.COORD_PRIMER && (
            <ExternalLink
              href="https://github.com/vector-engineering/covidcg/blob/master/static_data/primers.csv"
              style={{ marginLeft: '20px' }}
            >
              Primer/probe definitions
            </ExternalLink>
          )}
          {configStore.coordinateMode === COORDINATE_MODES.COORD_PRIMER && (
            <PrimerSelectContainer
              placeholderText={
                state.selectedPrimers.length === 0
                  ? 'Select or search...'
                  : state.selectedPrimers.length.toString() +
                    ' primers/probes selected...'
              }
            >
              {primerDropdown}
            </PrimerSelectContainer>
          )}
        </ModeRadioVertical>

        {/* CUSTOM COORDS */}
        <ModeRadioVertical>
          <ModeLabel>
            <input
              className="radio-input"
              type="radio"
              value={COORDINATE_MODES.COORD_CUSTOM}
              checked={
                configStore.coordinateMode === COORDINATE_MODES.COORD_CUSTOM
              }
              onChange={handleModeChange}
            />
            <span className="select-text">Custom Coordinates</span>
            {configStore.coordinateMode !== COORDINATE_MODES.COORD_CUSTOM && (
              <span className="hint-text">Select to show options</span>
            )}
            <UpdateButton
              show={
                state.customCoordinatesChanged &&
                configStore.coordinateMode === COORDINATE_MODES.COORD_CUSTOM
              }
              disabled={!state.validCustomCoords}
              onClick={handleCustomCoordSubmit}
            >
              Confirm
            </UpdateButton>
          </ModeLabel>
          {configStore.coordinateMode === COORDINATE_MODES.COORD_CUSTOM && (
            <CoordForm>
              <input
                type="text"
                value={state.customCoordText}
                onChange={handleCustomCoordChange}
              />
              <ReactTooltip
                className="filter-sidebar-tooltip"
                id="custom-coord-tooltip"
                type="light"
                effect="solid"
                border={true}
                borderColor="#888"
              />
              <QuestionButton
                data-tip='<p>Coordinates are in the form "start..end". Multiple ranges can be separated with ";"</p><p>i.e., "100..300;500..550"</p><p>Coordinates relative to the WIV04 reference sequence (EPI_ISL_402124)</p>'
                data-html="true"
                data-for="custom-coord-tooltip"
              />
            </CoordForm>
          )}
        </ModeRadioVertical>

        {/* CUSTOM SEQUENCES */}
        <ModeRadioVertical>
          <ModeLabel>
            <input
              className="radio-input"
              type="radio"
              value={COORDINATE_MODES.COORD_SEQUENCE}
              checked={
                configStore.coordinateMode === COORDINATE_MODES.COORD_SEQUENCE
              }
              onChange={handleModeChange}
            />
            <span className="select-text">Match Sequences</span>
            {configStore.coordinateMode !== COORDINATE_MODES.COORD_SEQUENCE && (
              <span className="hint-text">Select to show options</span>
            )}
            <UpdateButton
              show={
                state.customSequencesChanged &&
                configStore.coordinateMode === COORDINATE_MODES.COORD_SEQUENCE
              }
              disabled={!state.validCustomSequences}
              onClick={handleCustomSequencesSubmit}
            >
              Confirm
            </UpdateButton>
          </ModeLabel>
          {configStore.coordinateMode === COORDINATE_MODES.COORD_SEQUENCE && (
            <>
              <CoordForm>
                <ValidationInput
                  type="text"
                  value={state.customSequences}
                  onChange={handleCustomSequencesChange}
                  invalid={!state.validCustomSequences}
                />
                <ReactTooltip
                  className="filter-sidebar-tooltip"
                  id="custom-sequence-tooltip"
                  type="light"
                  effect="solid"
                  border={true}
                  borderColor="#888"
                />
                <QuestionButton
                  data-tip='<p>Select coordinates based on matches to the entered sequence (can be forward or reverse)</p><p>Please only enter A, T, C, or G. Enter in more than one sequence by separating them with ";"</p><p>Sequences are matched to the WIV04 reference sequence (EPI_ISL_402124)</p>'
                  data-html="true"
                  data-for="custom-sequence-tooltip"
                />
              </CoordForm>
              {!state.validCustomSequences && (
                <InvalidText>One or more sequences are invalid</InvalidText>
              )}
              {configStore.coordinateMode ===
                COORDINATE_MODES.COORD_SEQUENCE && (
                <RangesText>
                  Coordinates:{' '}
                  {configStore
                    .getCoordinateRanges()
                    .map((range) => range.join('..'))
                    .join(';')}
                </RangesText>
              )}
            </>
          )}
        </ModeRadioVertical>
      </ModeSelectForm>
    </SelectContainer>
  );
});

CoordinateSelect.propTypes = {};
CoordinateSelect.defaultProps = {};

export default CoordinateSelect;
