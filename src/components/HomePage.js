import React, { useState } from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import _ from 'underscore';

import GeneSelect from './GeneSelect';
import GroupBySelect from './GroupBySelect';
import DropdownContainer from './DropdownContainer';

import { VegaLite } from 'react-vega';

//import initial_entropy_spec from '../vega/barplot_v3.vl.json';
import areaStackSpecInitial from '../vega/area_stack.vl.json';

import { connect } from '../stores/connect';
import DataTable from './Table/DataTable';
import Header from './Header';
import SideBar from './Sidebar';
import { asyncStates } from '../stores/uiStore';
import SkeletonElement from './SkeletonElement';
import LoadingSpinner from './LoadingSpinner';

const HomePageDiv = styled.div`
  display: grid;
  grid-template-columns: [col1] 300px [col2] calc(100vw - 300px) [col3];
  grid-template-rows: [row1] 50px [row2] auto [row3];
  height: 100vh;
  width: 100vw;
`;
const FilterSidebar = styled.div`
  grid-column: col1 / col2;
  grid-row: row1 / row3;

  background-color: #f8f8f8;
  //padding-right: 10px;
  padding-bottom: 15px;
  border-right: 1px solid #aaa;
  display: flex;
  flex-direction: column;
`;
const PlotContainer = styled.div`
  grid-column: col2 / col3;
  grid-row: row2 / row3;

  width: 100%;
  box-sizing: border-box;

  padding-left: 10px;
  padding-top: 10px;

  .vega-embed {
    width: calc(100% - 110px);
  }
`;
const PlotOptions = styled.div``;

const AreaStackSelectContainer = styled.div`
  select {
    margin-left: 0.65em;
    padding: 1px 4px;
  }
`;

const AreaStackModeSelect = ({ mode, onChange }) => {
  return (
    <AreaStackSelectContainer>
      <label>
        Display mode:
        <select value={mode} onChange={onChange}>
          <option value="counts">Counts</option>
          <option value="percentages">Percentages</option>
        </select>
      </label>
    </AreaStackSelectContainer>
  );
};
AreaStackModeSelect.propTypes = {
  mode: PropTypes.string.isRequired,
  onChange: PropTypes.func.isRequired,
};

const HomePage = observer(({ covidStore, uiStore }) => {
  // 'percentages' or 'counts'
  const [areaStackMode, setAreaStackMode] = useState('percentages');

  const handleGroupingChange = (groupKey, dnaOrAa) => {
    covidStore.changeGrouping(groupKey, dnaOrAa);
  };

  const handleGeneChange = (gene) => {
    covidStore.selectGene(gene);
  };

  const handleBrush = (...args) => {
    //console.log(args);
    // this.setState({
    //   info: JSON.stringify(args),
    // });
    covidStore.selectDateRange(
      Object.prototype.hasOwnProperty.call(args[1], 'date')
        ? args[1].date
        : [-1, -1]
    );
  };

  const treeSelectOnChange = (currentNode, selectedNodes) => {
    covidStore.selectLocations(selectedNodes);
  };
  const treeSelectOnAction = (node, action) => {
    console.log('onAction::', action, node);
  };
  const treeSelectOnNodeToggleCurrentNode = (currentNode) => {
    console.log('onNodeToggle::', currentNode);
  };

  const onChangeAreaStackMode = (event) => setAreaStackMode(event.target.value);

  // Make a deep copy of the Vega spec so we can edit it
  const areaStackSpec = JSON.parse(JSON.stringify(areaStackSpecInitial));

  if (areaStackMode === 'percentages') {
    areaStackSpec['vconcat'][0]['encoding']['y']['stack'] = 'normalize';
  }

  // Adapt labels to groupings
  if (covidStore.groupKey === 'lineage') {
    // y-axis title
    areaStackSpec['vconcat'][0]['encoding']['y']['axis']['title'] =
      'Sequences by Lineage';
    // Tooltip title
    areaStackSpec['vconcat'][0]['encoding']['tooltip'][0]['title'] = 'Lineage';
  } else if (covidStore.groupKey === 'snp') {
    // y-axis title
    areaStackSpec['vconcat'][0]['encoding']['y']['axis']['title'] =
      'Sequences by ' + (covidStore.dnaOrAa === 'dna' ? 'NT' : 'AA') + ' SNP';
    // Tooltip title
    areaStackSpec['vconcat'][0]['encoding']['tooltip'][0]['title'] =
      (covidStore.dnaOrAa === 'dna' ? 'NT' : 'AA') + ' SNP';
  }

  const renderPlotContent = () => {
    if (uiStore.caseDataState === asyncStates.STARTED) {
      return (
        <div
          style={{
            paddingTop: '24px',
            paddingRight: '24px',
            paddingLeft: '12px',
          }}
        >
          <SkeletonElement delay={1} height={'436px'}>
            <LoadingSpinner />
          </SkeletonElement>
        </div>
      );
    } else {
      return (
        <VegaLite
          data={{
            case_data: covidStore.caseData,
          }}
          spec={areaStackSpec}
          signalListeners={{
            brush: _.debounce(handleBrush, 500),
          }}
        />
      );
    }
  };

  return (
    <>
      <HomePageDiv>
        <SideBar />
        <FilterSidebar>
          <Header />
          <GroupBySelect
            groupKey={covidStore.groupKey}
            dnaOrAa={covidStore.dnaOrAa}
            onChange={handleGroupingChange}
          />
          <GeneSelect
            genes={covidStore.genes}
            selectedGene={covidStore.selectedGene}
            onChange={handleGeneChange}
          />
          <DropdownContainer
            data={covidStore.selectTree.children}
            onChange={(currentNode, selectedNodes) => {
              treeSelectOnChange(currentNode, selectedNodes);
            }}
            onAction={treeSelectOnAction}
            onNodeToggle={treeSelectOnNodeToggleCurrentNode}
          />
        </FilterSidebar>

        <PlotContainer>
          {/* <VegaLite
          data={{
            entropy_data: entropy_data
          }} 
          spec={entropy_spec}
        /> */}
          <PlotOptions>
            <AreaStackModeSelect
              mode={areaStackMode}
              onChange={onChangeAreaStackMode}
            />
          </PlotOptions>
          {renderPlotContent()}
          <DataTable />
        </PlotContainer>
      </HomePageDiv>
    </>
  );
});

HomePage.propTypes = {
  covidStore: PropTypes.object.isRequired,
  router: PropTypes.object.isRequired,
};

// eslint-disable-next-line react/display-name
export default connect(HomePage);
