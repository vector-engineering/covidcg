import React from 'react';
// import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import ReactTooltip from 'react-tooltip';
import QuestionButton from '../Buttons/QuestionButton';

import { SelectContainer, RadioForm, Link } from './GroupBySelect.styles';

import {
  GROUP_SNV,
  DNA_OR_AA,
  COORDINATE_MODES,
} from '../../constants/defs.json';
import { config } from '../../config';

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
  } else if (
    configStore.coordinateMode === COORDINATE_MODES.COORD_GENE &&
    configStore.selectedGene.name === 'All Genes'
  ) {
    aaDisabled = true;
    aaDisabledMessage = ' (please select one gene)';
  } else if (
    configStore.coordinateMode === COORDINATE_MODES.COORD_PROTEIN &&
    configStore.selectedProtein.name === 'All Proteins'
  ) {
    aaDisabled = true;
    aaDisabledMessage = ' (please select one protein)';
  } else if (
    configStore.coordinateMode === COORDINATE_MODES.COORD_GENE &&
    configStore.selectedGene.protein_coding === 0
  ) {
    aaDisabled = true;
    aaDisabledMessage = ' (please select protein-coding gene)';
  }

  const groupSelectItems = [];
  Object.keys(config.group_cols).forEach((group) => {
    groupSelectItems.push(
      <div className="radio-item" key={`select-${group}`}>
        <label>
          <input
            className="radio-input"
            type="radio"
            value={group}
            checked={configStore.groupKey === group}
            onChange={handleGroupKeyChange}
          />
          <span>{config.group_cols[group].title}</span>
        </label>
      </div>
    );
  });

  return (
    <SelectContainer>
      <ReactTooltip
        className="filter-sidebar-tooltip"
        id="groupby-tooltip"
        type="light"
        effect="solid"
        border={true}
        borderColor="#888"
      />
      <RadioForm>
        <span className="form-title">
          Group sequences by
          <QuestionButton
            data-tip={`<p>For "SNV", count SNVs <i>independently</i> across all selected sequences. i.e., the sum of all SNV counts may exceed the total number of selected sequences, as sequences may have multiple SNVs.</p><p>For any other mode, e.g., "Lineage", aggregate sequences by their "Lineage" assignment. In this case, the sum of all lineage counts will equal the total number of selected sequences, as a sequence has one and only one "Lineage" assignment`}
            data-html="true"
            data-for="groupby-tooltip"
          />
        </span>
        <div className="radio-row">
          {groupSelectItems}
          <div className="radio-item">
            <label>
              <input
                className="radio-input"
                type="radio"
                value={GROUP_SNV}
                checked={configStore.groupKey === GROUP_SNV}
                onChange={handleGroupKeyChange}
              />
              <span>SNV</span>
            </label>
          </div>
        </div>
        {Object.keys(config.group_cols).includes(configStore.groupKey) && (
          <>
            {config.group_cols[configStore.groupKey].description}
            <Link href={config.group_cols[configStore.groupKey].link.href}>
              {config.group_cols[configStore.groupKey].link.title}
            </Link>
          </>
        )}
      </RadioForm>
      <RadioForm>
        <span className="form-title">Mutation format</span>
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