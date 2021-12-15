import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';

import QuestionButton from '../Buttons/QuestionButton';

import {
  SelectContainer,
  RadioForm,
  Link,
  HintText,
} from './GroupBySelect.styles';

import {
  GROUP_MUTATION,
  DNA_OR_AA,
  COORDINATE_MODES,
} from '../../constants/defs.json';
import { config } from '../../config';

const GroupBySelect = observer(
  ({
    groupKey,
    dnaOrAa,
    coordinateMode,
    selectedGene,
    selectedProtein,
    selectedReference,
    onGroupKeyChange,
    onDnaOrAaChange,
    onReferenceChange,
  }) => {
    let handleGroupKeyChange = (event) => {
      onGroupKeyChange(event.target.value);
    };

    let handleDnaOrAaChange = (event) => {
      onDnaOrAaChange(event.target.value);
    };

    let handleReferenceChange = (event) => {
      onReferenceChange(event.target.value);
    };

    let aaDisabledMessage = '';
    let aaDisabled = false;
    if (
      coordinateMode !== COORDINATE_MODES.COORD_GENE &&
      coordinateMode !== COORDINATE_MODES.COORD_PROTEIN
    ) {
      aaDisabledMessage = ' (only for gene/protein)';
      aaDisabled = true;
    } else if (
      coordinateMode === COORDINATE_MODES.COORD_GENE &&
      selectedGene.name === 'All Genes'
    ) {
      aaDisabled = true;
      aaDisabledMessage = ' (please select one gene)';
    } else if (
      coordinateMode === COORDINATE_MODES.COORD_PROTEIN &&
      selectedProtein.name === 'All Proteins'
    ) {
      aaDisabled = true;
      aaDisabledMessage = ' (please select one protein)';
    } else if (
      coordinateMode === COORDINATE_MODES.COORD_GENE &&
      selectedGene.protein_coding === 0
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
              checked={groupKey === group}
              onChange={handleGroupKeyChange}
            />
            <span>{config.group_cols[group].title}</span>
          </label>
        </div>
      );
    });

    const renderGroupKeySelect = () => {
      return (
        <>
          <div className="radio-row">
            {groupSelectItems}
            <div className="radio-item">
              <label>
                <input
                  className="radio-input"
                  type="radio"
                  value={GROUP_MUTATION}
                  checked={groupKey === GROUP_MUTATION}
                  onChange={handleGroupKeyChange}
                />
                <span>Mutation</span>
              </label>
            </div>
          </div>
          {Object.keys(config.group_cols).includes(groupKey) && (
            <>
              {config.group_cols[groupKey].description}
              <Link href={config.group_cols[groupKey].link.href}>
                {config.group_cols[groupKey].link.title}
              </Link>
            </>
          )}
        </>
      );
    };

    const renderDnaOrAaSelect = () => {
      return (
        <div className="radio-row">
          <div className="radio-item">
            <input
              type="radio"
              id="dnaChoice"
              name="dnaOrAa"
              value={DNA_OR_AA.DNA}
              checked={dnaOrAa === DNA_OR_AA.DNA}
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
              checked={dnaOrAa === DNA_OR_AA.AA}
              disabled={aaDisabled}
              onChange={handleDnaOrAaChange}
            ></input>
            <label htmlFor="aaChoice">AA</label>
            <span className="disabled-text">{aaDisabledMessage}</span>
          </div>
        </div>
      );
    };

    const renderRSVRefSelect = () => {
      return (
        <div className="radio-row">
          <div className="radio-item">
            <input
              type="radio"
              id="RSVAChoice"
              name="rsvAorB"
              value="A"
              checked={selectedReference === 'A'}
              onChange={handleReferenceChange}
            ></input>
            <label htmlFor="RSVAChoice">RSV-A</label>
          </div>
          <div className="radio-item">
            <input
              type="radio"
              id="RSVBChoice"
              name="rsvAorB"
              value="B"
              checked={selectedReference === 'B'}
              onChange={handleReferenceChange}
            ></input>
            <label htmlFor="RSVAChoice">RSV-B</label>
          </div>
        </div>
      );
    };

    return (
      <SelectContainer>
        <RadioForm>
          <span className="form-title">
            Group sequences by
            <QuestionButton
              data-tip={`<p>For "Mutation", count mutations <i>independently</i> across all selected sequences. i.e., the sum of all mutation counts may exceed the total number of selected sequences, as sequences may have multiple mutations.</p><p>For any other mode, e.g., "Lineage", aggregate sequences by their "Lineage" assignment. In this case, the sum of all lineage counts will equal the total number of selected sequences, as a sequence has one and only one "Lineage" assignment</p>`}
              data-html="true"
              data-place="right"
              data-for="main-tooltip"
            />
          </span>
          {renderGroupKeySelect()}
        </RadioForm>
        <RadioForm>
          <span className="form-title">Mutation format</span>
          {groupKey !== GROUP_MUTATION && (
            <HintText>
              Switch to &quot;Mutation&quot; under &quot;Group sequences
              by&quot; in order to enable Mutation Formatting
            </HintText>
          )}
          {groupKey === GROUP_MUTATION && renderDnaOrAaSelect()}
        </RadioForm>
        {config.virus === 'rsv' && (
          <RadioForm>
            <span className="form-title">Reference Sequence</span>
            {renderRSVRefSelect()}
          </RadioForm>
        )}
      </SelectContainer>
    );
  }
);

GroupBySelect.propTypes = {
  // groupKey: PropTypes.string,
  // dnaOrAa: PropTypes.string,
  // coordinateMode: PropTypes.string,
  // selectedGene: PropTypes.object,
  // selectedProtein: PropTypes.object
  onGroupKeyChange: PropTypes.func,
  onDnaOrAaChange: PropTypes.func,
  onReferenceChange: PropTypes.func,
};
GroupBySelect.defaultProps = {};

export default GroupBySelect;
