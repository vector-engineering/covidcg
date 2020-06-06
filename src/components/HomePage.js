import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';

import _ from 'underscore';

import GeneSelect from './GeneSelect';
import DropdownContainer from './DropdownContainer';

import 'react-dropdown-tree-select/dist/styles.css';

import { VegaLite } from 'react-vega';

//import initial_entropy_spec from '../vega/barplot_v3.vl.json';
import area_stack_absolute_spec from '../vega/area_stack.vl.json';
import area_stack_norm_spec from '../vega/area_stack_norm.vl.json';

import '../styles/home-page.scss';
import { connect } from '../stores/connect';
import LineageDataTable from './LineageDataTable';
import Header from './Header';

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

    this.handleGeneChange = this.handleGeneChange.bind(this);

    this.treeSelectOnChange = this.treeSelectOnChange.bind(this);
    this.treeSelectOnAction = this.treeSelectOnAction.bind(this);
    this.treeSelectOnNodeToggleCurrentNode = this.treeSelectOnNodeToggleCurrentNode.bind(
      this
    );

    this.onChangeAreaStackMode = this.onChangeAreaStackMode.bind(this);
  }

  handleGeneChange(gene) {
    console.log('Gene change:', gene);

    this.props.covidStore.selectGene(event.target.value);
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
      <div className="home-page">
        <div className="filter-sidebar">
          <GeneSelect
            genes={covidStore.genes}
            value={covidStore.selectedGene}
            startPos={covidStore.startPos}
            endPos={covidStore.endPos}
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
        </div>
        <Header />
        <div className="plot-container">
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
        </div>
      </div>
    );
  }
}

HomePage.propTypes = {
  covidStore: PropTypes.object.isRequired,
  router: PropTypes.object.isRequired,
};

// eslint-disable-next-line react/display-name
export default connect(HomePage);
