import React from 'react';
import PropTypes from 'prop-types';
import { NavLink } from 'react-router-dom';
import _ from 'underscore';
import { connect } from 'react-redux';
import { bindActionCreators } from 'redux';
import * as actions from '../actions/covidActions';

import { getReferenceSequence } from '../utils/lineageData';

import GeneSelect from './GeneSelect';
import DropdownContainer from './DropdownContainer';
import HeatmapCell from './HeatmapCell';
import DataTable from 'react-data-table-component';

import 'react-dropdown-tree-select/dist/styles.css';

import { VegaLite } from 'react-vega';

//import initial_entropy_spec from '../vega/barplot_v3.vl.json';
import area_stack_absolute_spec from '../vega/area_stack.vl.json';
import area_stack_norm_spec from '../vega/area_stack_norm.vl.json';

import '../styles/home-page.scss';
import LiteMolCell from './Cells/LiteMolCell';

export class HomePage extends React.Component {
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

    this.props.actions.selectGene(event.target.value);
  }

  handleBrush(...args) {
    //console.log(args);
    // this.setState({
    //   info: JSON.stringify(args),
    // });
    this.props.actions.selectDateRange(
      Object.prototype.hasOwnProperty.call(args[1], 'date')
        ? args[1].date
        : [-1, -1]
    );
  }

  treeSelectOnChange(currentNode, selectedNodes) {
    console.log('onChange::', currentNode, selectedNodes);
    this.props.actions.selectLocations(selectedNodes);
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
    const activeStyle = { color: 'blue' };

    // Get the bases at the positions, for the reference sequence
    let ref_seq = getReferenceSequence();

    // Build a column for each changing position
    let pos_cols = [];
    this.props.covid.changingPositions.forEach((pos) => {
      pos_cols.push({
        name: pos.toString(),
        selector: 'pos_' + pos.toString(),
        sortable: false,
        width: '40px',
        center: true,
        compact: true,
        style: {
          fontFamily: 'monospace',
          fontWeight: '500',
          fontSize: '1.25em',
        },
        conditionalCellStyles: [
          {
            when: (row) => row['pos_' + pos.toString()] != ref_seq[pos],
            style: {
              backgroundColor: '#FFFF00',
            },
          },
        ],
      });
    });

    let maxCasesSum = _.reduce(
      this.props.covid.caseDataAggLineageList,
      (memo, lineage) => Math.max(memo, lineage.cases_sum),
      0
    );
    let minCasesSum = _.reduce(
      this.props.covid.caseDataAggLineageList,
      (memo, lineage) => Math.min(memo, lineage.cases_sum),
      0
    );
    let maxCasesPercent = _.reduce(
      this.props.covid.caseDataAggLineageList,
      (memo, lineage) => Math.max(memo, lineage.cases_percent),
      0
    );
    let minCasesPercent = _.reduce(
      this.props.covid.caseDataAggLineageList,
      (memo, lineage) => Math.min(memo, lineage.cases_percent),
      0
    );

    //console.log(this.props.covid);

    let area_stack_spec =
      this.state.area_stack_mode === 'percentages'
        ? area_stack_norm_spec
        : area_stack_absolute_spec;

    return (
      <div className="home-page">
        <div className="filter-sidebar">
          <GeneSelect
            genes={this.props.covid.genes}
            value={this.props.covid.selectedGene}
            startPos={this.props.covid.startPos}
            endPos={this.props.covid.endPos}
            onChange={this.handleGeneChange}
          />
          <DropdownContainer
            data={this.props.covid.selectTree.children}
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
        <div className="header">
          <div className="title-container">
            <h1>COVID-UI</h1>
          </div>
          <div className="nav-links">
            <NavLink exact to="/" activeStyle={activeStyle}>
              Home
            </NavLink>
            <NavLink to="/about" activeStyle={activeStyle}>
              About
            </NavLink>
            <a
              href="https://github.com/vector-engineering/covid_ui"
              target="_blank"
              rel="noopener noreferrer"
            >
              View on GitHub
            </a>
          </div>
        </div>
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
              case_data: this.props.covid.caseData,
            }}
            spec={area_stack_spec}
            signalListeners={this.handlers}
          />

          <DataTable
            className="data-table"
            data={this.props.covid.caseDataAggLineageList}
            columns={[
              {
                name: 'Lineage',
                selector: 'lineage',
                sortable: true,
                width: '100px',
                style: {
                  fontWeight: '700',
                },
              },
              {
                name: 'Cases',
                selector: 'cases_sum',
                sortable: true,
                width: '85px',
                cell: (row) => {
                  return (
                    <HeatmapCell
                      value={row.cases_sum}
                      min={minCasesSum}
                      max={maxCasesSum}
                      percent={false}
                    />
                  );
                },
              },
              {
                name: '% Cases',
                selector: 'cases_percent',
                sortable: true,
                width: '85px',
                cell: (row) => {
                  return (
                    <HeatmapCell
                      value={row.cases_percent}
                      min={minCasesPercent}
                      max={maxCasesPercent}
                      percent={true}
                    />
                  );
                },
              },
              //   {
              //     name: 'jmol',
              //     selector: 'jmol',
              //     sortable: false,
              //     width: '100%',
              //     cell: () => {
              //       return <LiteMolCell />;
              //     },
              //   },
            ].concat(pos_cols)}
            striped={true}
            highlightOnHover={true}
            dense={true}
            // fixedHeader={true}
            // fixedHeaderScrollHeight={'400px'}

            pagination={false}
            defaultSortField={'lineage'}
            defaultSortAsc={true}
            conditionalRowStyles={[
              {
                when: (row) => row.lineage == 'Reference',
                style: 'background-color: #dff3fe !important;',
              },
            ]}
            sortFunction={(rows, field, direction) => {
              // Set aside the reference, and remove it from the rows list
              let refRow = _.findWhere(rows, { lineage: 'Reference' });
              rows = _.reject(rows, (row) => row.lineage == 'Reference');

              // Normal sorting...
              rows = _.sortBy(rows, (row) => {
                return row[field];
              });
              // Reverse if descending
              if (direction == 'desc') {
                rows.reverse();
              }
              // Add the reference row to the beginning
              rows.unshift(refRow);

              return rows;
            }}
            customStyles={{
              headCells: {
                style: {
                  paddingLeft: '8px', // override the cell padding for head cells
                  paddingRight: '8px',
                },
              },
              cells: {
                style: {
                  paddingLeft: '8px', // override the cell padding for data cells
                  paddingRight: '8px',
                },
              },
            }}
          />
        </div>
      </div>
    );
  }
}

HomePage.propTypes = {
  actions: PropTypes.object.isRequired,
  covid: PropTypes.object.isRequired,
};

function mapStateToProps(state) {
  return {
    covid: state.covid,
  };
}

function mapDispatchToProps(dispatch) {
  return {
    actions: bindActionCreators(actions, dispatch),
  };
}

export default connect(mapStateToProps, mapDispatchToProps)(HomePage);
