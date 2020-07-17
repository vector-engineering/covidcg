import React from 'react';
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

const ModeSelectForm = styled.form`
  .radio {
    margin-left: 10px;
    margin-bottom: 5px;

    label.checkbox-label {
      display: flex;
      flex-direction: row;
      align-items: center;

      padding-right: 10px;

      input.checkbox-input {
        margin: 0px 8px 0px 0px;
      }
    }
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
const genes = getAllGenes();
const proteins = getAllProteins();

const GeneSelect = observer(() => {
  const { covidStore } = useStores();

  const changeCoordinateMode = ({
    coordinateMode,
    selectedGene,
    selectedProtein,
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
        <div className="radio">
          <label className="checkbox-label">
            <input
              className="checkbox-input"
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
          </label>
        </div>
        <div className="radio">
          <label className="checkbox-label">
            <input
              className="checkbox-input"
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
          </label>
        </div>
      </ModeSelectForm>
    </SelectContainer>
  );
});

GeneSelect.propTypes = {};
GeneSelect.defaultProps = {};

export default GeneSelect;
