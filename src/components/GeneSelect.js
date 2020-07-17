import React, { useState, useEffect } from 'react';
// import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../stores/connect';
import styled from 'styled-components';

import { getAllGenes } from '../utils/gene';
import { getAllProteins } from '../utils/protein';

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

const ModeSelectForm = styled.form``;

const ModeRadio = styled.div`
  margin-left: 10px;
  margin-bottom: 5px;
`;

const ModeHorizontal = styled.label`
  display: flex;
  flex-direction: row;
  align-items: center;

  padding-right: 10px;

  input.radio-input {
    margin: 0px 8px 0px 0px;
  }
`;

const ModeVertical = styled.label`
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
    margin-left: 0.65em;
    padding: 1px 5px;
    width: 100%;
    border-radius: 3px;
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

const GeneSelect = observer(() => {
  const { covidStore } = useStores();

  const [state, setState] = useState({
    customStart: covidStore.customCoordinates[0],
    customEnd: covidStore.customCoordinates[1],
  });

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
    changeCoordinateMode({ selectedGene: event.target.value });
  };

  const handleProteinChange = (event) => {
    changeCoordinateMode({ selectedProtein: event.target.value });
  };

  const handleCustomCoordStartChange = (event) => {
    setState({ ...state, customStart: event.target.value });
  };
  const handleCustomCoordEndChange = (event) => {
    setState({ ...state, customEnd: event.target.value });
  };
  const handleCustomCoordSubmit = (event) => {
    event.preventDefault();
    changeCoordinateMode({
      customCoordinates: [state.customStart, state.customEnd],
    });
  };

  // Create option elements

  // GENE
  let geneOptionElements = [];
  // All Genes option
  geneOptionElements.push(
    <option key="all" value="all">
      All Genes
    </option>
  );

  genes.forEach((gene) => {
    geneOptionElements.push(
      <option key={gene.gene} value={gene.gene}>
        {gene.gene}&nbsp;&nbsp;({gene.start}..{gene.end})
      </option>
    );
  });

  // PROTEIN
  let proteinOptionElements = [];
  // All Proteins option
  proteinOptionElements.push(
    <option key="all" value="all">
      All Proteins
    </option>
  );

  proteins.forEach((protein) => {
    proteinOptionElements.push(
      <option key={protein.protein} value={protein.protein}>
        {protein.protein}&nbsp;&nbsp;({protein.segments})
      </option>
    );
  });

  return (
    <SelectContainer>
      <span className="title">Genomic Coordinates</span>
      <ModeSelectForm>
        <ModeRadio>
          <ModeHorizontal>
            <input
              className="radio-input"
              type="radio"
              value="gene"
              checked={covidStore.coordinateMode === 'gene'}
              onChange={handleModeChange}
            />
            <SelectForm>
              <label>Gene:</label>
              <select
                value={covidStore.selectedGene.gene}
                onChange={handleGeneChange}
              >
                {geneOptionElements}
              </select>
            </SelectForm>
          </ModeHorizontal>
        </ModeRadio>
        <ModeRadio>
          <ModeHorizontal>
            <input
              className="radio-input"
              type="radio"
              value="protein"
              checked={covidStore.coordinateMode === 'protein'}
              onChange={handleModeChange}
            />
            <SelectForm>
              <label>Protein:</label>
              <select
                value={covidStore.selectedProtein.protein}
                onChange={handleProteinChange}
              >
                {proteinOptionElements}
              </select>
            </SelectForm>
          </ModeHorizontal>
        </ModeRadio>
        <ModeRadio>
          <ModeVertical>
            <input
              className="radio-input"
              type="radio"
              value="custom"
              checked={covidStore.coordinateMode === 'custom'}
              onChange={handleModeChange}
            />
            Custom Coordinates:
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
          </ModeVertical>
        </ModeRadio>
      </ModeSelectForm>
    </SelectContainer>
  );
});

GeneSelect.propTypes = {};
GeneSelect.defaultProps = {};

export default GeneSelect;
