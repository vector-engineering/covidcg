import React, { useState, useEffect, useMemo } from 'react';
// import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../stores/connect';
import styled from 'styled-components';
import _ from 'underscore';

import Button from './Buttons/Button';
import DropdownTreeSelect from 'react-dropdown-tree-select';

import { getAllGenes } from '../utils/gene';
import { getAllProteins } from '../utils/protein';
import {
  getPrimerSelectTree,
  getPrimerByName,
  getPrimersByGroup,
} from '../utils/primer';

const SelectContainer = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;
  justify-content: flex-start;

  margin: 5px 5px 5px 5px;
  padding: 0px 8px 5px 8px;

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
    margin-right: 5px;
    max-width: 4rem;
  }
  button {
    flex-grow: 1;
  }
`;

const genes = getAllGenes();
const proteins = getAllProteins();

const CoordinateSelect = observer(() => {
  const { covidStore } = useStores();

  // Create option elements

  // GENE
  let geneOptionElements = [];
  genes.forEach((gene) => {
    geneOptionElements.push(
      <option key={gene.gene} value={gene.gene}>
        {gene.gene}&nbsp;&nbsp;({gene.start}..{gene.end})
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
    customStart: covidStore.customCoordinates[0],
    customEnd: covidStore.customCoordinates[1],
  });

  // Disable "All Genes" and "All Proteins" option
  // when in AA mode and non-SNP grouping
  // useEffect(() => {
  //   let _geneOptionElements = state.geneOptionElements;
  //   let _proteinOptionElements = state.proteinOptionElements;

  //   if (covidStore.groupKey !== 'snp' && covidStore.dnaOrAa === 'aa') {

  //   }

  // }, [covidStore.groupKey, covidStore.dnaOrAa]);

  // Update custom coordinates from the store
  useEffect(() => {
    setState({
      ...state,
      customStart: covidStore.customCoordinates[0],
      customEnd: covidStore.customCoordinates[1],
    });
  }, [covidStore.customCoordinates]);

  const changeCoordinateMode = ({
    coordinateMode,
    selectedGene,
    selectedProtein,
    selectedPrimers,
    customCoordinates,
  }) => {
    covidStore.changeCoordinateMode({
      coordinateMode:
        coordinateMode === undefined
          ? covidStore.coordinateMode
          : coordinateMode,
      selectedGene:
        selectedGene === undefined
          ? covidStore.selectedGene.gene
          : selectedGene,
      selectedProtein:
        selectedProtein === undefined
          ? covidStore.selectedProtein.protein
          : selectedProtein,
      selectedPrimers:
        selectedPrimers === undefined
          ? covidStore.selectedPrimers
          : selectedPrimers,
      customCoordinates:
        customCoordinates === undefined
          ? covidStore.customCoordinates
          : customCoordinates,
    });
  };

  const handleModeChange = (event) => {
    changeCoordinateMode({ coordinateMode: event.target.value });
  };

  const handleGeneChange = (event) => {
    changeCoordinateMode({
      coordinateMode: 'gene',
      selectedGene: event.target.value,
    });
  };

  const handleProteinChange = (event) => {
    changeCoordinateMode({
      coordinateMode: 'protein',
      selectedProtein: event.target.value,
    });
  };

  const handleCustomCoordStartChange = (event) => {
    setState({ ...state, customStart: event.target.value });
  };
  const handleCustomCoordEndChange = (event) => {
    setState({ ...state, customEnd: event.target.value });
  };
  const handleCustomCoordSubmit = (event) => {
    event.preventDefault();
    // Change to custom mode implicitly
    changeCoordinateMode({
      coordinateMode: 'custom',
      customCoordinates: [state.customStart, state.customEnd],
    });
  };

  const checkPrimersChanged = (selectedPrimers) => {
    // Is the current covidStore.selectedPrimers the same as the current selection?
    let changed = selectedPrimers.length !== covidStore.selectedPrimers.length;
    if (!changed && selectedPrimers.length > 0) {
      // Run through once - both selectedPrimers and the covidStore version are sorted
      for (let i = 0; i < selectedPrimers.length; i++) {
        if (!_.isEqual(selectedPrimers[i], covidStore.selectedPrimers[i])) {
          changed = true;
          break;
        }
      }
    }
    return changed;
  };

  useEffect(() => {
    setState({
      ...state,
      primersChanged: checkPrimersChanged(state.selectedPrimers),
    });
  }, [covidStore.selectedPrimers]);

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
      coordinateMode: 'primer',
      selectedPrimers: state.selectedPrimers,
    });
  };

  // This component needs to be in a memoized function
  // since it manages its own local state. It should never be re-rendered
  // forcefully
  const primerDropdown = useMemo(
    () => (
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
    ),
    [state.primerTreeData]
  );

  return (
    <SelectContainer>
      <span className="title">Genomic Coordinates</span>
      <ModeSelectForm>
        <ModeRadioHorizontal>
          <ModeLabel>
            <input
              className="radio-input"
              type="radio"
              value="gene"
              checked={covidStore.coordinateMode === 'gene'}
              onChange={handleModeChange}
            />
            Gene
          </ModeLabel>
          <SelectForm>
            <select
              value={covidStore.selectedGene.gene}
              onChange={handleGeneChange}
            >
              <option
                key="All Genes"
                value="All Genes"
                disabled={
                  covidStore.groupKey !== 'snp' && covidStore.dnaOrAa === 'aa'
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
              value="protein"
              checked={covidStore.coordinateMode === 'protein'}
              onChange={handleModeChange}
            />
            Protein
          </ModeLabel>
          <SelectForm>
            <select
              value={covidStore.selectedProtein.protein}
              onChange={handleProteinChange}
            >
              <option
                key="All Proteins"
                value="All Proteins"
                disabled={
                  covidStore.groupKey !== 'snp' && covidStore.dnaOrAa === 'aa'
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
              value="primer"
              checked={covidStore.coordinateMode === 'primer'}
              onChange={handleModeChange}
            />
            Primers/Probes:
          </ModeLabel>
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
              value="custom"
              checked={covidStore.coordinateMode === 'custom'}
              onChange={handleModeChange}
            />
            Custom Coordinates:
          </ModeLabel>
          <CoordForm>
            <span>From</span>
            <input
              type="number"
              min={1}
              max={29903}
              step={1}
              value={state.customStart}
              onChange={handleCustomCoordStartChange}
            />
            <span>To</span>
            <input
              type="number"
              min={1}
              max={29903}
              step={1}
              value={state.customEnd}
              onChange={handleCustomCoordEndChange}
            />
            <button onClick={handleCustomCoordSubmit}>OK</button>
          </CoordForm>
        </ModeRadioVertical>
      </ModeSelectForm>
    </SelectContainer>
  );
});

CoordinateSelect.propTypes = {};
CoordinateSelect.defaultProps = {};

export default CoordinateSelect;
