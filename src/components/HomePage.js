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
import LineageDataTable from './LineageDataTable';
import Header from './Header';
import SideBar from './Sidebar';

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
  box-shadow: 0px 0px 5px #aaa;
  display: block;

  .location-tree-title {
    margin-left: 12px;
  }
`;
const PlotContainer = styled.div`
  grid-column: col2 / col3;
  grid-row: row2 / row3;

  width: 100%;
  box-sizing: border-box;

  padding-left: 10px;
  padding-top: 10px;

  .vega-embed {
    width: calc(100% - 100px);
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

const HomePage = observer(({ covidStore }) => {
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
    console.log('onChange::', currentNode, selectedNodes);
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

  // Adapt axis labels to groupings
  if (covidStore.groupKey === 'lineage') {
    areaStackSpec['vconcat'][0]['encoding']['y']['axis']['title'] =
      'Cases by Lineage';
  } else if (covidStore.groupKey === 'snp') {
    areaStackSpec['vconcat'][0]['encoding']['y']['axis']['title'] =
      'Cases by ' + (covidStore.dnaOrAa === 'dna' ? 'NT' : 'AA') + ' SNP';
  }

  console.log(areaStackSpec);

  return (
    <HomePageDiv>
      <SideBar />
      <FilterSidebar>
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
        <span className="location-tree-title">Selected Locations:</span>
        <DropdownContainer
          data={covidStore.selectTree.children}
          onChange={treeSelectOnChange}
          onAction={treeSelectOnAction}
          onNodeToggle={treeSelectOnNodeToggleCurrentNode}
        />
      </FilterSidebar>
      <Header />
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

        <VegaLite
          data={{
            case_data: covidStore.caseData,
          }}
          spec={areaStackSpec}
          signalListeners={{
            brush: _.debounce(handleBrush, 500),
          }}
        />

        <LineageDataTable />
      </PlotContainer>
    </HomePageDiv>
  );
});

HomePage.propTypes = {
  covidStore: PropTypes.object.isRequired,
  router: PropTypes.object.isRequired,
};

// eslint-disable-next-line react/display-name
export default connect(HomePage);
