import React, { useState, useEffect, useMemo } from 'react';
import { observer } from 'mobx-react';

import ExternalLink from '../Common/ExternalLink';
import QuestionButton from '../Buttons/QuestionButton';
import {
  PrimerSelect,
  ModeLabel,
  SelectContainer,
  ModeSelectForm,
  ModeRadioVertical,
  SelectForm,
  CoordForm,
  PrimerSelectContainer,
  ValidationInput,
  InvalidText,
  RangesText,
  DomainSelectForm,
  HintText,
} from './CoordinateSelect.styles';

import { coordsToText, textToCoords } from '../../utils/coordinates';
import { getAllGenes, getAllProteins } from '../../utils/gene_protein';
import {
  getPrimerSelectTree,
  getPrimerByName,
  getPrimersByGroup,
} from '../../utils/primer';
import { queryReferenceSequence, getReference } from '../../utils/reference';
import { config } from '../../config';

import {
  DNA_OR_AA,
  COORDINATE_MODES,
  GROUP_MUTATION,
} from '../../constants/defs.json';

const CoordinateSelect = observer(
  ({
    groupKey,
    dnaOrAa,
    coordinateMode,
    selectedGene,
    selectedProtein,
    selectedPrimers,
    selectedReference,
    residueCoordinates,
    validResidueCoordinates,
    customCoordinates,
    validCustomCoordinates,
    customSequences,
    validCustomSequences,

    updateCoordinateMode,
    updateSelectedGene,
    updateSelectedProtein,
    updateResidueCoordinates,
    updateValidResidueCoordinates,
    updateSelectedPrimers,
    updateCustomCoordinates,
    updateValidCustomCoordinates,
    updateCustomSequences,
    updateValidCustomSequences,
  }) => {
    // Create option elements
    let genes = {};
    let proteins = {};

    genes = getAllGenes(selectedReference);
    proteins = getAllProteins(selectedReference);

    // GENE
    let geneOptionElements = [];
    genes.forEach((gene) => {
      geneOptionElements.push(
        <option
          key={gene.name}
          value={gene.name}
          disabled={!gene.protein_coding && dnaOrAa === DNA_OR_AA.AA}
        >
          {gene.name}&nbsp;&nbsp;(
          {gene.segments.map((range) => range.join('..')).join(';')})
        </option>
      );
    });

    // GENE DOMAINS
    let geneDomainOptionElements = {};
    genes.forEach((gene) => {
      if (!gene.protein_coding) {
        return;
      }
      geneDomainOptionElements[gene.name] = [
        <option
          key={`${gene.name}-default`}
          value={`${gene.name}-default`}
          disabled={true}
        >
          - select an option -
        </option>,
        <option key={`${gene.name}-all`} value={`${gene.name}-all`}>
          Entire {gene.name} Gene ({gene.residue_offset_range.join('..')})
        </option>,
      ];
      gene.domains.forEach((domain, j) => {
        geneDomainOptionElements[gene.name].push(
          <option key={`${gene.name}-${domain.name}-${j}`} value={domain.name}>
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
        <option key={protein.name} value={protein.name}>
          {protein.name}&nbsp;&nbsp;(
          {protein.segments.map((range) => range.join('..')).join(';')})
        </option>
      );
    });

    // PROTEIN DOMAINS
    let proteinDomainOptionElements = {};
    proteins.forEach((protein) => {
      proteinDomainOptionElements[protein.name] = [
        <option
          key={`${protein.name}-default`}
          value={`${protein.name}-default`}
          disabled={true}
        >
          {' '}
          - select an option -
        </option>,
        <option key={`${protein.name}-all`} value={`${protein.name}-all`}>
          Entire {protein.name} Protein (
          {protein.residue_offset_range.join('..')})
        </option>,
      ];
      protein.domains.forEach((domain, j) => {
        proteinDomainOptionElements[protein.name].push(
          <option
            key={`${protein.name}-${domain.name}-${j}`}
            value={domain.name}
          >
            {domain.name}&nbsp;&nbsp;(
            {domain.ranges.map((range) => range.join('..')).join(';')})
          </option>
        );
      });
    });

    const [state, setState] = useState({
      primerTreeData: Object.assign(getPrimerSelectTree()),

      customCoordText: coordsToText(customCoordinates),
      customSequences: customSequences.join(';'),

      residueCoordsText: residueCoordinates
        .map((range) => range.join('..'))
        .join(';'),
    });

    // Disable "All Genes" and "All Proteins" option
    // when in AA mode and non-mutation grouping
    // useEffect(() => {
    //   let _geneOptionElements = state.geneOptionElements;
    //   let _proteinOptionElements = state.proteinOptionElements;
    //   if (groupKey !== GROUP_MUTATION && dnaOrAa === DNA_OR_AA.AA) {
    //   }
    // }, [groupKey, dnaOrAa]);

    // Update custom coordinates from the store
    useEffect(() => {
      setState({
        ...state,
        customCoordText: coordsToText(customCoordinates),
      });
    }, [customCoordinates]);

    // Update custom sequences from the store
    useEffect(() => {
      setState({
        ...state,
        customSequences: customSequences.join(';'),
      });
    }, [customSequences]);

    const handleModeChange = (event) => {
      updateCoordinateMode(event.target.value);
    };

    const handleGeneChange = (event) => {
      updateSelectedGene(event.target.value);
    };

    const handleProteinChange = (event) => {
      updateSelectedProtein(event.target.value);
    };

    // Use a regex to match numbers, since just because JS
    // can parse an integer, doesn't mean it should...
    const numPattern = /^([0-9-]+)$/;

    const handleResidueCoordsChange = (event) => {
      // Parse current custom coordinates
      const curResidueCoords = event.target.value
        .split(';')
        .map((range) => range.split('..'));

      // Check that these are valid
      const validResidueCoordinates = !curResidueCoords.some((range) => {
        // Return true if invalid
        return (
          range.length !== 2 ||
          numPattern.exec(range[0]) === null ||
          numPattern.exec(range[1]) === null ||
          parseInt(range[0]) > parseInt(range[1]) ||
          parseInt(range[0]) < 1 - selectedGene.residue_offset ||
          parseInt(range[1]) > selectedGene.len_aa - selectedGene.residue_offset
        );
      });

      setState({
        ...state,
        residueCoordsText: event.target.value,
      });

      if (validResidueCoordinates) {
        updateResidueCoordinates(
          curResidueCoords.map((range) => range.map((coord) => parseInt(coord)))
        );
      } else {
        updateValidResidueCoordinates(validResidueCoordinates);
      }
    };

    // Use the selected domain to fill in the residue coordinates input
    const handleGeneDomainChange = (event) => {
      const domainName = event.target.value;
      const newResidueCoords = [];
      let newResidueCoordsText;

      if (event.target.value === selectedGene.name + '-all') {
        newResidueCoordsText = `${selectedGene.residue_offset_range[0]}..${selectedGene.residue_offset_range[1]}`;
        newResidueCoords.push(selectedGene.residue_offset_range);
      } else {
        const domainObj = selectedGene.domains.find(
          (domain) => domain.name === domainName
        );

        domainObj.ranges.forEach((range) =>
          newResidueCoords.push(range.slice())
        );
        newResidueCoordsText = domainObj.ranges
          .map((range) => range.join('..'))
          .join(';');
      }

      setState({
        ...state,
        residueCoordsText: newResidueCoordsText,
      });

      updateResidueCoordinates(newResidueCoords);
    };

    const handleProteinDomainChange = (event) => {
      const domainName = event.target.value;
      const newResidueCoords = [];
      let newResidueCoordsText;

      if (event.target.value === selectedProtein.name + '-all') {
        newResidueCoordsText = `${selectedProtein.residue_offset_range[0]}..${selectedProtein.residue_offset_range[1]}`;
        newResidueCoords.push(selectedProtein.residue_offset_range);
      } else {
        const domainObj = selectedProtein.domains.find(
          (domain) => domain.name === domainName
        );

        domainObj.ranges.forEach((range) =>
          newResidueCoords.push(range.slice())
        );

        newResidueCoordsText = domainObj.ranges
          .map((range) => range.join('..'))
          .join(';');
      }

      setState({
        ...state,
        residueCoordsText: newResidueCoordsText,
      });

      updateResidueCoordinates(newResidueCoords);
    };

    // Update residue coordinates from store
    useEffect(() => {
      setState({
        ...state,
        residueCoordsText: residueCoordinates
          .map((range) => range.join('..'))
          .join(';'),
      });
    }, [residueCoordinates]);

    const handleCustomCoordChange = (event) => {
      // Parse current custom coordinates
      const curCustomCoords = textToCoords(event.target.value);
      // Check that these are valid
      const validCustomCoordinates =
        !curCustomCoords.some((range) => {
          // Return true if invalid
          return (
            range.length !== 3 || // Each range must consist of 3 elements
            !config.segments.includes(range[0]) || // Segment must be valid
            numPattern.exec(range[1]) === null || // Start/end must be integers
            numPattern.exec(range[2]) === null ||
            parseInt(range[1]) > parseInt(range[2]) || // Start cannot be greater than end
            // Range must be within [1, segment sequence length]
            range[1] < 1 ||
            range[2] >
              getReference(selectedReference).segments[range[0]]['sequence']
                .length
          );
        }) && new Set(curCustomCoords.map((range) => range[0])).size === 1; // Segment must be the same for all ranges
      setState({
        ...state,
        customCoordText: event.target.value,
      });
      if (validCustomCoordinates) {
        updateCustomCoordinates(
          curCustomCoords.map((range) => {
            range[1] = parseInt(range[1]);
            range[2] = parseInt(range[2]);
            return range;
          })
        );
      } else {
        updateValidCustomCoordinates(validCustomCoordinates);
      }
    };

    const handleCustomSequencesChange = (event) => {
      const curText = event.target.value.toUpperCase();
      const sequences = curText.split(';');
      // Check that the query sequence is not empty, and that the
      // reference sequence includes the sequence,
      const validCustomSequences = !sequences.some((seq) => {
        return (
          seq.length === 0 ||
          queryReferenceSequence(selectedReference, seq) === 0
        );
      });

      setState({
        ...state,
        customSequences: curText,
      });

      if (validCustomSequences) {
        updateCustomSequences(curText.split(';'));
      } else {
        updateValidCustomSequences(validCustomSequences);
      }
    };

    // Check all selected primers
    useEffect(() => {
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

      selectedPrimers.forEach((primer) => {
        const institutionNode = primerTreeData.find(
          (node) => node.value === primer.Institution
        );
        if (institutionNode === undefined) {
          return;
        }
        const primerNode = institutionNode.children.find(
          (child) => child.value === primer.Name
        );
        if (primerNode === undefined) {
          return;
        }
        primerNode.checked = true;
      });

      setState({
        ...state,
        primerTreeData,
      });
    }, [selectedPrimers]);

    const onPrimerSelect = (currentNode, selectedNodes) => {
      let selectedPrimers = [];
      selectedNodes.forEach((node) => {
        if (node.level === 'group') {
          selectedPrimers = selectedPrimers.concat(
            getPrimersByGroup(node.value)
          );
        } else if (node.level === 'individual') {
          selectedPrimers = selectedPrimers.concat(getPrimerByName(node.value));
        }
      });

      // Sort by Institution then Name
      selectedPrimers = selectedPrimers
        .map((primer) => {
          primer.sortKey = primer.Institution.concat('-', primer.Name);
          return primer;
        })
        .sort((a, b) => {
          return a.sortKey > b.sortKey;
        });

      updateSelectedPrimers(selectedPrimers);
    };

    // Maintain tree expansion state
    const onPrimerTreeNodeToggle = (currentNode) => {
      const primerTreeData = state.primerTreeData.slice();

      primerTreeData.forEach((node) => {
        if (node.value === currentNode.value) {
          node.expanded = currentNode.expanded;
        }
      });

      setState({
        ...state,
        primerTreeData,
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
          showDropdown="always"
          texts={{
            placeholder: 'Search...',
            noMatches: 'No matches found',
          }}
          onChange={onPrimerSelect}
          onNodeToggle={onPrimerTreeNodeToggle}
        />
      );
    }, [state.primerTreeData]);

    const renderMainForm = () => {
      return (
        <ModeSelectForm>
          {/* GENE SELECT */}
          <ModeRadioVertical>
            <ModeLabel>
              <input
                className="radio-input"
                type="radio"
                value={COORDINATE_MODES.COORD_GENE}
                checked={coordinateMode === COORDINATE_MODES.COORD_GENE}
                onChange={handleModeChange}
              />
              <span className="option-text">Gene</span>
              <SelectForm>
                <select value={selectedGene.name} onChange={handleGeneChange}>
                  <option
                    key="All Genes"
                    value="All Genes"
                    disabled={dnaOrAa === DNA_OR_AA.AA}
                  >
                    All Genes
                  </option>
                  {geneOptionElements}
                </select>
              </SelectForm>
            </ModeLabel>
            {coordinateMode === COORDINATE_MODES.COORD_GENE &&
              selectedGene.name !== 'All Genes' &&
              selectedGene.protein_coding && (
                <>
                  <CoordForm>
                    <span className="coord-prefix">Residue indices:</span>
                    <input
                      type="text"
                      value={state.residueCoordsText}
                      onChange={handleResidueCoordsChange}
                    />
                    <QuestionButton
                      rebuildAfterMount={true}
                      data-tip='
                        <p>
                          Coordinates are in the form "start..end". 
                          Multiple ranges can be separated with ";"
                        </p>
                        <p>
                          i.e., "100..300;500..550"
                        </p>
                        <p>
                          Coordinates are relative to the gene ORF
                        </p>'
                      data-html="true"
                      data-for="main-tooltip"
                    />
                  </CoordForm>
                  {!validResidueCoordinates && (
                    <InvalidText>Invalid coordinate format</InvalidText>
                  )}
                  <DomainSelectForm>
                    <span>Domain:</span>
                    <select
                      value={`${selectedGene.name}-default`}
                      onChange={handleGeneDomainChange}
                    >
                      {geneDomainOptionElements[selectedGene.name]}
                    </select>
                    <QuestionButton
                      rebuildAfterMount={true}
                      data-tip='
                        <p>
                          Coordinates relative to the gene ORF, and are in the form "start..end".
                        </p>
                        <p>
                          Selecting a domain will replace the range(s) to the residue indices input
                        </p>'
                      data-html="true"
                      data-for="main-tooltip"
                    />
                  </DomainSelectForm>
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
                checked={coordinateMode === COORDINATE_MODES.COORD_PROTEIN}
                onChange={handleModeChange}
              />
              <span className="option-text">Protein</span>
              <SelectForm>
                <select
                  value={selectedProtein.name}
                  onChange={handleProteinChange}
                >
                  <option
                    key="All Proteins"
                    value="All Proteins"
                    disabled={dnaOrAa === DNA_OR_AA.AA}
                  >
                    All Proteins
                  </option>
                  {proteinOptionElements}
                </select>
              </SelectForm>
            </ModeLabel>
            {coordinateMode === COORDINATE_MODES.COORD_PROTEIN &&
              selectedProtein.name !== 'All Proteins' && (
                <>
                  <CoordForm>
                    <span className="coord-prefix">Residue indices:</span>
                    <input
                      type="text"
                      value={state.residueCoordsText}
                      onChange={handleResidueCoordsChange}
                    />
                    <QuestionButton
                      rebuildAfterMount={true}
                      data-tip='
                        <p>
                          Coordinates are in the form "start..end". 
                          Multiple ranges can be separated with ";"
                        </p>
                        <p>
                          i.e., "100..300;500..550"
                        </p>
                        <p>
                          Coordinates are relative to the protein ORF
                        </p>'
                      data-html="true"
                      data-for="main-tooltip"
                    />
                  </CoordForm>
                  {!validResidueCoordinates && (
                    <InvalidText>Invalid coordinate format</InvalidText>
                  )}
                  <DomainSelectForm>
                    <span>Domain:</span>
                    <select
                      value={`${selectedProtein.name}-default`}
                      onChange={handleProteinDomainChange}
                    >
                      {proteinDomainOptionElements[selectedProtein.name]}
                    </select>
                    <QuestionButton
                      rebuildAfterMount={true}
                      data-tip='
                        <p>
                          Coordinates relative to the protein ORF, and are in 
                          the form "start..end".
                        </p>
                        <p>
                          Selecting a domain will replace the range(s) to the residue indices input
                        </p>'
                      data-html="true"
                      data-for="main-tooltip"
                    />
                  </DomainSelectForm>
                </>
              )}
          </ModeRadioVertical>

          {/* PRIMER/PROBE SELECT */}
          {config.virus === 'sars2' && (
            <ModeRadioVertical>
              <ModeLabel>
                <input
                  className="radio-input"
                  type="radio"
                  value={COORDINATE_MODES.COORD_PRIMER}
                  checked={coordinateMode === COORDINATE_MODES.COORD_PRIMER}
                  onChange={handleModeChange}
                />
                <span className="select-text">Primers/Probes</span>
                {coordinateMode !== COORDINATE_MODES.COORD_PRIMER && (
                  <span className="hint-text">Select to show options</span>
                )}
              </ModeLabel>
              {coordinateMode === COORDINATE_MODES.COORD_PRIMER && (
                <ExternalLink
                  href="https://github.com/vector-engineering/covidcg/blob/master/static_data/sars2/primers.csv"
                  style={{ marginLeft: '20px' }}
                >
                  Primer/probe definitions
                </ExternalLink>
              )}
              {coordinateMode === COORDINATE_MODES.COORD_PRIMER && (
                <PrimerSelectContainer
                  placeholderText={
                    selectedPrimers.length === 0
                      ? 'Select or search...'
                      : selectedPrimers.length.toString() +
                        ' primers/probes selected...'
                  }
                >
                  {primerDropdown}
                </PrimerSelectContainer>
              )}
            </ModeRadioVertical>
          )}

          {/* CUSTOM COORDS */}
          <ModeRadioVertical>
            <ModeLabel>
              <input
                className="radio-input"
                type="radio"
                value={COORDINATE_MODES.COORD_CUSTOM}
                checked={coordinateMode === COORDINATE_MODES.COORD_CUSTOM}
                onChange={handleModeChange}
              />
              <span className="select-text">Custom Coordinates</span>
              {coordinateMode !== COORDINATE_MODES.COORD_CUSTOM && (
                <span className="hint-text">Select to show options</span>
              )}
            </ModeLabel>
            {coordinateMode === COORDINATE_MODES.COORD_CUSTOM && (
              <CoordForm>
                <input
                  type="text"
                  value={state.customCoordText}
                  onChange={handleCustomCoordChange}
                />
                <QuestionButton
                  rebuildAfterMount={true}
                  // Show segment formatting only if this viruses has more than one segment
                  data-tip={`
                    <p>
                      Coordinates are in the form ${
                        config.segments.length > 1
                          ? '"segment:start..end"'
                          : '"start..end"'
                      }. 
                      Multiple ranges can be separated with ";".
                    </p>
                    ${
                      config.segments.length > 1
                        ? '<p>If no segment is defined, it will default to the first segment. Note: currently, only one segment at a time can be specified. i.e., you cannot specify ranges from or across multiple segments</p>'
                        : ''
                    }
                    ${
                      config.segments.length > 1
                        ? '<p>i.e., "2:100..300;2:500..550"</p>'
                        : '<p>i.e., "100..300;500..550"</p>'
                    }
                    <p>
                      Coordinates relative to the reference genome: <b>${selectedReference}</b>
                    </p>`}
                  data-html="true"
                  data-for="main-tooltip"
                />
              </CoordForm>
            )}
            {!validCustomCoordinates && (
              <InvalidText>Invalid coordinate format</InvalidText>
            )}
          </ModeRadioVertical>

          {/* CUSTOM SEQUENCES */}
          <ModeRadioVertical>
            <ModeLabel>
              <input
                className="radio-input"
                type="radio"
                value={COORDINATE_MODES.COORD_SEQUENCE}
                checked={coordinateMode === COORDINATE_MODES.COORD_SEQUENCE}
                onChange={handleModeChange}
              />
              <span className="select-text">Match Sequences</span>
              {coordinateMode !== COORDINATE_MODES.COORD_SEQUENCE && (
                <span className="hint-text">Select to show options</span>
              )}
            </ModeLabel>
            {coordinateMode === COORDINATE_MODES.COORD_SEQUENCE && (
              <>
                <CoordForm>
                  <ValidationInput
                    type="text"
                    value={state.customSequences}
                    onChange={handleCustomSequencesChange}
                    invalid={!validCustomSequences}
                  />
                  <QuestionButton
                    rebuildAfterMount={true}
                    data-tip={`
                      <p>
                        Select coordinates based on matches to the entered sequence 
                        (can be forward or reverse)
                      </p>
                      <p>
                        Please only enter A, T, C, or G. 
                        Enter in more than one sequence by separating them with ";"
                      </p>
                      <p>
                        Sequences are matched to the reference genome: <b>${selectedReference}</b>
                      </p>`}
                    data-html="true"
                    data-for="main-tooltip"
                  />
                </CoordForm>
                {!validCustomSequences && (
                  <InvalidText>One or more sequences are invalid</InvalidText>
                )}
                {validCustomSequences &&
                  coordinateMode === COORDINATE_MODES.COORD_SEQUENCE && (
                    <RangesText>
                      Coordinates:{' '}
                      {coordsToText(
                        state.customSequences
                          .split(';')
                          .map((seq) =>
                            queryReferenceSequence(selectedReference, seq)
                          )
                      )}
                    </RangesText>
                  )}
              </>
            )}
          </ModeRadioVertical>
        </ModeSelectForm>
      );
    };

    return (
      <SelectContainer>
        <span className="title">
          Genomic Coordinates
          <QuestionButton
            data-tip='
              <p>
                When grouping by mutation, only show mutations within the given genomic coordinates.
              </p>
              <p>
                When grouping by lineage/clade, only show consensus mutations 
                within the given genomic coordinates.
              </p>
              <p>
                These options are only enabled when in "Mutation" mode.
              </p>'
            data-html="true"
            data-for="main-tooltip"
          />
        </span>
        {groupKey !== GROUP_MUTATION && (
          <ModeSelectForm>
            <HintText>
              Switch to &quot;Mutation&quot; under &quot;Group sequences
              by&quot; in order to enable Genomic Coordinate filtering.
            </HintText>
          </ModeSelectForm>
        )}
        {groupKey === GROUP_MUTATION && renderMainForm()}
      </SelectContainer>
    );
  }
);

CoordinateSelect.propTypes = {};
CoordinateSelect.defaultProps = {};

export default CoordinateSelect;
