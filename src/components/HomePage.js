import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import _ from 'underscore';

import GeneSelect from './GeneSelect';
import GroupBySelect from './GroupBySelect';
import DropdownContainer from './DropdownContainer';

import { VegaLite } from 'react-vega';

//import initial_entropy_spec from '../vega/barplot_v3.vl.json';
import area_stack_absolute_spec from '../vega/area_stack.vl.json';
import area_stack_norm_spec from '../vega/area_stack_norm.vl.json';

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
`;
const PlotContainer = styled.div`
  grid-column: col2 / col3;
  grid-row: row2 / row3;

  padding-left: 10px;
  padding-top: 10px;
`;

@observer
class HomePage extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {
      area_stack_mode: 'percentages',
    };
    //this.processEntropyData.bind(this);

    this.handleBrush = this.handleBrush.bind(this);
    this.handlers = { brush: _.debounce(this.handleBrush, 500) };

    this.handleGroupingChange = this.handleGroupingChange.bind(this);
    this.handleGeneChange = this.handleGeneChange.bind(this);

    this.treeSelectOnChange = this.treeSelectOnChange.bind(this);
    this.treeSelectOnAction = this.treeSelectOnAction.bind(this);
    this.treeSelectOnNodeToggleCurrentNode = this.treeSelectOnNodeToggleCurrentNode.bind(
      this
    );

    this.onChangeAreaStackMode = this.onChangeAreaStackMode.bind(this);
  }

  handleGroupingChange(groupKey, dnaOrAa) {
    this.props.covidStore.changeGrouping(groupKey, dnaOrAa);
  }

  handleGeneChange(gene) {
    this.props.covidStore.selectGene(gene);
  }

  handleBrush(...args) {
    const { covidStore } = this.props;
    //console.log(args);
    // this.setState({
    //   info: JSON.stringify(args),
    // });
    covidStore.selectDateRange(
      Object.prototype.hasOwnProperty.call(args[1], 'date')
        ? args[1].date
        : [-1, -1]
    );
  }

  treeSelectOnChange(currentNode, selectedNodes) {
    console.log('onChange::', currentNode, selectedNodes);
    this.props.covidStore.selectLocations(selectedNodes);
  }
  treeSelectOnAction(node, action) {
    console.log('onAction::', action, node);
  }
  treeSelectOnNodeToggleCurrentNode(currentNode) {
    console.log('onNodeToggle::', currentNode);
  }

  onChangeAreaStackMode(e) {
    this.setState({
      area_stack_mode: e.target.value,
    });
  }

  render() {
    const { covidStore } = this.props;
    const activeStyle = { color: 'blue' };

    //console.log(.covid);

    let area_stack_spec =
      this.state.area_stack_mode === 'percentages'
        ? area_stack_norm_spec
        : area_stack_absolute_spec;

    return (
      <HomePageDiv>
        <SideBar />
        <FilterSidebar>
          <GroupBySelect
            groupKey={covidStore.groupKey}
            dnaOrAa={covidStore.dnaOrAa}
            onChange={this.handleGroupingChange}
          />
          <GeneSelect
            genes={covidStore.genes}
            selectedGene={covidStore.selectedGene}
            onChange={this.handleGeneChange}
          />
          <DropdownContainer
            data={covidStore.selectTree.children}
            onChange={this.treeSelectOnChange}
            onAction={this.treeSelectOnAction}
            onNodeToggle={this.treeSelectOnNodeToggleCurrentNode}
            className="geo-dropdown-tree-select"
            clearSearchOnChange={false}
            keepTreeOnSearch={true}
            keepChildrenOnSearch={true}
            showPartiallySelected={true}
            showDropdown="always"
            inlineSearchInput={true}
            texts={{
              placeholder: 'Choose...',
              noMatches: 'No matches found',
            }}
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

          <div className="plot-options">
            <label>
              Display mode
              <select
                value={this.state.area_stack_mode}
                onChange={this.onChangeAreaStackMode}
              >
                <option value="counts">Counts</option>
                <option value="percentages">Percentages</option>
              </select>
            </label>
          </div>

          <VegaLite
            data={{
              case_data: covidStore.caseData,
            }}
            spec={area_stack_spec}
            signalListeners={this.handlers}
          />

          <LineageDataTable />
        </PlotContainer>
      </HomePageDiv>
    );
  }
}

HomePage.propTypes = {
  covidStore: PropTypes.object.isRequired,
  router: PropTypes.object.isRequired,
};

// eslint-disable-next-line react/display-name
export default connect(HomePage);
