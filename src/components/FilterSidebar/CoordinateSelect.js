import React, { useState, useEffect, useMemo } from 'react';
// import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import styled from 'styled-components';
import _ from 'underscore';

import ExternalLink from '../Common/ExternalLink';
import Button from '../Buttons/Button';
import DropdownTreeSelect from 'react-dropdown-tree-select';
import QuestionButton from '../Buttons/QuestionButton';

import { getAllGenes } from '../../utils/gene';
import { getAllProteins } from '../../utils/protein';
import {
  getPrimerSelectTree,
  getPrimerByName,
  getPrimersByGroup,
} from '../../utils/primer';
import { referenceSequenceIncludes } from '../../utils/reference';

import {
  GROUP_KEYS,
  DNA_OR_AA,
  COORDINATE_MODES,
} from '../../constants/config';

const SelectContainer = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;
  justify-content: flex-start;

  margin: 5px 5px 0px 5px;
  padding: 0px 8px 0px 8px;

  span.title {
    margin-bottom: 5px;
  }
`;

const ModeSelectForm = styled.div``;

const ModeRadioHorizontal = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;

  margin-left: 10px;
  margin-bottom: 5px;
`;

const ModeRadioVertical = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;
  justify-content: flex-start;

  margin-left: 10px;
  margin-bottom: 8px;
`;

const ModeLabel = styled.label`
  display: flex;
  flex-direction: row;
  align-items: center;

  padding-right: 10px;

  input.radio-input {
    margin: 0px 8px 0px 0px;
  }
`;

const SelectForm = styled.form`
  display: flex;
  flex-direction: row;
  align-items: center;
  flex-grow: 1;

  label {
    display: flex;
    flex-direction: row;
    align-items: center;
    justify-content: flex-start;
  }

  select {
    background-color: white;
    flex-grow: 1;
    padding: 1px 5px;
    width: 100%;
    border-radius: 3px;
  }
`;

const UpdatePrimersButton = styled(Button)`
  display: ${(props) => (props.show ? 'block' : 'none')};
  font-size: 1em;
  margin-left: 20px;
`;

UpdatePrimersButton.defaultProps = {
  show: false,
};

const PrimerSelectContainer = styled.div`
  span.placeholder {
    &:after {
      content: "${(props) => props.placeholderText}";
    }
  }
`;

PrimerSelectContainer.defaultProps = {
  placeholderText: '',
};

const PrimerSelect = styled(DropdownTreeSelect)`
  a.dropdown-trigger {
    &:focus {
      outline: none;
    }

    ul.tag-list {
      margin-top: 6px;
      margin-bottom: 0px;
      list-style: none;
      padding-left: 20px;

      li.tag-item {
        margin-right: 10px;
        line-height: normal;

        span.placeholder {
          background-color: #ffffff;
          font-size: 0em;
          &:after {
            font-size: 0.9rem;
            font-weight: normal;
            border: 1px solid #888;
            border-radius: 3px;
            background-color: #fff;
            padding: 3px 8px;
          }
          &:hover,
          &:focus {
            &:after {
              background-color: #eee;
              border: 1px solid #666;
            }
          }
        }
        span.tag {
          display: none;
        }
      }
    }
  }

  .dropdown-content {
    display: flex;
    flex-direction: column;
    align-items: stretch;

    padding: 5px;
    padding-top: 8px;

    input.search {
      padding: 3px 8px;
      margin-left: 15px;
      margin-right: 5px;
      border: 1px solid #888;
      border-radius: 3px;
      font-size: 1em;
      &:focus {
        outline: none;
      }
    }

    ul.root {
      margin: 5px 0px;
      padding-left: 15px;
      font-weight: normal;
      font-size: 1em;

      li.node {
        i.toggle {
          font-family: monospace;
          font-size: 1.25em;
          font-style: normal;
          font-weight: 500;
          &:hover {
            color: #888888;
          }

          white-space: pre;
          margin-right: 4px;
          outline: none;

          cursor:pointer &:after {
            content: ' ';
          }
          &.collapsed:after {
            content: '+';
          }
          &.expanded:after {
            content: '-';
          }
        }

        label {
          .checkbox-item,
          .radio-item {
            vertical-align: middle;
            margin: 0 4px 0 0;
            &.simple-select {
              display: none;
            }
          }
        }
      }
    }
  }
`;

const CoordForm = styled.form`
  display: flex;
  flex-direction: row;
  align-items: center;

  margin-top: 3px;
  margin-left: 20px;
  padding-right: 10px;

  span {
    font-weight: normal;
    margin-right: 5px;
  }
  input {
    flex-grow: 1;
    margin-right: 5px;
    max-width: 15rem;
  }
  button {
    flex-grow: 1;
  }
`;

const UpdateButton = styled(Button)`
  display: ${(props) => (props.show ? 'block' : 'none')};
  font-size: 0.9em;
  padding: 3px 8px;
  margin-left: 10px;
`;
UpdateButton.defaultProps = {
  show: false,
};

const ValidationInput = styled.input`
  border: 1px solid ${({ invalid }) => (invalid ? '#dc3545' : '#aaa')};
  &:focus {
    border: 2px solid ${({ invalid }) => (invalid ? '#dc3545' : '#aaa')};
  }
`;
ValidationInput.defaultProps = {
  invalid: false,
};

const InvalidText = styled.span`
  margin-left: 20px;
  font-size: 0.9em;
  font-weight: normal;
  line-height: normal;
  color: #dc3545;
`;

const RangesText = styled.span`
  margin-left: 20px;
  font-size: 0.9em;
  font-weight: normal;
  line-height: normal;
  color: #888;
`;

const genes = getAllGenes();
const proteins = getAllProteins();

const CoordinateSelect = observer(() => {
  const { configStore } = useStores();

  // Create option elements

  // GENE
  let geneOptionElements = [];
  genes.forEach((gene) => {
    geneOptionElements.push(
      <option key={gene.gene} value={gene.gene}>
        {gene.gene}&nbsp;&nbsp;({gene.segments})
      </option>
    );
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
      validCustomCoords,
      customCoordinatesChanged: storeCustomCoords !== event.target.value,
      customCoordText: event.target.value,
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
        <ModeRadioHorizontal>
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
            Gene
          </ModeLabel>
          <SelectForm>
            <select
              value={configStore.selectedGene.gene}
              onChange={handleGeneChange}
            >
              <option
                key="All Genes"
                value="All Genes"
                disabled={
                  configStore.groupKey !== GROUP_KEYS.GROUP_SNV &&
                  configStore.dnaOrAa === DNA_OR_AA.AA
                }
              >
                All Genes
              </option>
              {geneOptionElements}
            </select>
          </SelectForm>
        </ModeRadioHorizontal>
        <ModeRadioHorizontal>
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
            Protein
          </ModeLabel>
          <SelectForm>
            <select
              value={configStore.selectedProtein.protein}
              onChange={handleProteinChange}
            >
              <option
                key="All Proteins"
                value="All Proteins"
                disabled={
                  configStore.groupKey !== GROUP_KEYS.GROUP_SNV &&
                  configStore.dnaOrAa === DNA_OR_AA.AA
                }
              >
                All Proteins
              </option>
              {proteinOptionElements}
            </select>
          </SelectForm>
        </ModeRadioHorizontal>

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
            <span>Primers/Probes</span>
          </ModeLabel>
          {configStore.coordinateMode === COORDINATE_MODES.COORD_PRIMER && (
            <ExternalLink
              href="https://github.com/vector-engineering/covidcg/blob/master/static_data/primers.csv"
              style={{ marginLeft: '20px' }}
            >
              (Primer/probe definitions)
            </ExternalLink>
          )}
          <UpdatePrimersButton
            show={state.primersChanged}
            onClick={updatePrimerSelection}
          >
            Update Primer Selection
          </UpdatePrimersButton>
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
        </ModeRadioVertical>

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
            <span>Custom Coordinates</span>
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
          <CoordForm>
            <input
              type="text"
              value={state.customCoordText}
              onChange={handleCustomCoordChange}
            />
            <QuestionButton
              data-tip='<p>Coordinates are in the form "start..end". Multiple ranges can be separated with ";"</p><p>i.e., "100..300;500..550"</p><p>Coordinates relative to Wuhan-Hu-1 reference sequence (NC_045512.2)</p>'
              data-html="true"
              data-for="tooltip-filter-sidebar"
            />
          </CoordForm>
        </ModeRadioVertical>
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
            <span>Match Sequences</span>
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
          <CoordForm>
            <ValidationInput
              type="text"
              value={state.customSequences}
              onChange={handleCustomSequencesChange}
              invalid={!state.validCustomSequences}
            />
            <QuestionButton
              data-tip='<p>Select coordinates based on matches to the entered sequence (can be forward or reverse)</p><p>Please only enter A, T, C, or G. Enter in more than one sequence by separating them with ";"</p><p>Sequences are matched to Wuhan-Hu-1 reference sequence (NC_045512.2)</p>'
              data-html="true"
              data-for="tooltip-filter-sidebar"
            />
          </CoordForm>
          {!state.validCustomSequences && (
            <InvalidText>One or more sequences are invalid</InvalidText>
          )}
          {configStore.coordinateMode === COORDINATE_MODES.COORD_SEQUENCE && (
            <RangesText>
              Coordinates:{' '}
              {configStore
                .getCoordinateRanges()
                .map((range) => range.join('..'))
                .join(';')}
            </RangesText>
          )}
        </ModeRadioVertical>
      </ModeSelectForm>
    </SelectContainer>
  );
});

CoordinateSelect.propTypes = {};
CoordinateSelect.defaultProps = {};

export default CoordinateSelect;
