import React from 'react';
// import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import styled from 'styled-components';

import ExternalLink from '../Common/ExternalLink';

import {
  GROUP_KEYS,
  DNA_OR_AA,
  COORDINATE_MODES,
} from '../../constants/config';

const SelectContainer = styled.div`
  padding: 7px 13px 0px 13px;
  box-shadow: 0px 2px 4px #aaa;
  background-color: #fcfcfc;
`;

const RadioForm = styled.form`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  margin-bottom: 10px;
  // margin-bottom: 5px;

  // background-color: #fff;
  // padding: 5px 10px;
  // border: 1px solid #ccc;

  span {
    flex-shrink: 0;
    margin-right: 10px;
  }

  .radio-row {
    display: flex;
    flex-direction: row;
    align-items: center;
    justify-content: flex-start;

    margin-top: 3px;
    margin-left: 12px;
  }

  .radio-item {
    margin-left: 0.65em;
    &:first-child {
      margin-left: 0px;
    }

    display: flex;
    flex-direction: row;
    align-items: center;

    input {
      margin: 0px;
      margin-right: 0.5em;
    }

    span.disabled-text {
      font-weight: normal;
      font-size: 0.8em;
      color: #888;
      flex-shrink: 1;
      white-space: normal;
      line-height: normal;
      margin-left: 5px;
    }
  }
`;

const Link = styled(ExternalLink)`
  font-size: 0.9em;
  margin-left: 10px;
`;
const CladeText = styled.span`
  font-size: 0.9em;
  margin-left: 10px;
  line-height: normal;
  font-weight: normal;
`;

const GroupBySelect = observer(() => {
  const { configStore } = useStores();

  let handleGroupKeyChange = (event) => {
    configStore.changeGrouping(event.target.value, configStore.dnaOrAa);
  };

  let handleDnaOrAaChange = (event) => {
    configStore.changeGrouping(configStore.groupKey, event.target.value);
  };

  let aaDisabledMessage = '';
  let aaDisabled = false;
  if (
    configStore.coordinateMode !== COORDINATE_MODES.COORD_GENE &&
    configStore.coordinateMode !== COORDINATE_MODES.COORD_PROTEIN
  ) {
    aaDisabledMessage = ' (only for gene/protein)';
    aaDisabled = true;
  } else if (configStore.groupKey !== GROUP_KEYS.GROUP_SNV) {
    if (
      configStore.coordinateMode === COORDINATE_MODES.COORD_GENE &&
      configStore.selectedGene.gene === 'All Genes'
    ) {
      aaDisabled = true;
      aaDisabledMessage = ' (please select one gene)';
    } else if (
      configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN &&
      configStore.selectedProtein.protein === 'All Proteins'
    ) {
      aaDisabled = true;
      aaDisabledMessage = ' (please select one protein)';
    }
  }

  const renderLineageLink = () => {
    if (configStore.groupKey === GROUP_KEYS.GROUP_LINEAGE) {
      return (
        <Link href="https://cov-lineages.org/descriptions.html">
          (Lineage Descriptions)
        </Link>
      );
    }

    return null;
  };

  const renderCladeLink = () => {
    if (configStore.groupKey === GROUP_KEYS.GROUP_CLADE) {
      return (
        <CladeText>
          For more information about clade and lineage nomenclature, visit this{' '}
          <ExternalLink href="https://www.gisaid.org/references/statements-clarifications/clade-and-lineage-nomenclature-aids-in-genomic-epidemiology-of-active-hcov-19-viruses/">
            [GISAID note]
          </ExternalLink>
        </CladeText>
      );
    }

    return null;
  };

  return (
    <SelectContainer>
      <RadioForm>
        <span className="form-title">Group sequences by</span>
        <div className="radio-row">
          <div className="radio-item">
            <label>
              <input
                className="radio-input"
                type="radio"
                value={GROUP_KEYS.GROUP_LINEAGE}
                checked={configStore.groupKey === GROUP_KEYS.GROUP_LINEAGE}
                onChange={handleGroupKeyChange}
              />
              <span>Lineage</span>
            </label>
          </div>
          <div className="radio-item">
            <label>
              <input
                className="radio-input"
                type="radio"
                value={GROUP_KEYS.GROUP_CLADE}
                checked={configStore.groupKey === GROUP_KEYS.GROUP_CLADE}
                onChange={handleGroupKeyChange}
              />
              <span>Clade</span>
            </label>
          </div>
          <div className="radio-item">
            <label>
              <input
                className="radio-input"
                type="radio"
                value={GROUP_KEYS.GROUP_SNV}
                checked={configStore.groupKey === GROUP_KEYS.GROUP_SNV}
                onChange={handleGroupKeyChange}
              />
              <span>SNV</span>
            </label>
          </div>
        </div>
        {renderLineageLink()}
        {renderCladeLink()}
      </RadioForm>
      <RadioForm>
        <span>Mutation format</span>
        <div className="radio-row">
          <div className="radio-item">
            <input
              type="radio"
              id="dnaChoice"
              name="dnaOrAa"
              value={DNA_OR_AA.DNA}
              checked={configStore.dnaOrAa === DNA_OR_AA.DNA}
              onChange={handleDnaOrAaChange}
            ></input>
            <label htmlFor="dnaChoice">NT</label>
          </div>
          <div className="radio-item">
            <input
              type="radio"
              id="aaChoice"
              name="dnaOrAa"
              value={DNA_OR_AA.AA}
              checked={configStore.dnaOrAa === DNA_OR_AA.AA}
              disabled={aaDisabled}
              onChange={handleDnaOrAaChange}
            ></input>
            <label htmlFor="aaChoice">AA</label>
            <span className="disabled-text">{aaDisabledMessage}</span>
          </div>
        </div>
      </RadioForm>
    </SelectContainer>
  );
});

GroupBySelect.propTypes = {};
GroupBySelect.defaultProps = {};

export default GroupBySelect;
