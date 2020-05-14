import React, { Component } from 'react';
import PropTypes from 'prop-types'
import DropdownTreeSelect from 'react-dropdown-tree-select';
import _ from 'underscore';

class DropdownContainer extends Component {
  constructor(props) {
    super(props);
    this.state = { 
        data: props.data 
    };
  }

  UNSAFE_componentWillReceiveProps = (nextProps) => {
    if(!_.isEqual(nextProps.data, this.state.data)) {
      this.setState({ data: nextProps.data });
    }
  }

  shouldComponentUpdate = (nextProps) => {
    return !_.isEqual(nextProps.data, this.state.data);
  }

  render() {
    const { data, ...rest } = this.props;
    data;
    return (
      <DropdownTreeSelect data={this.state.data} {...rest} />
    )
  }
}

DropdownContainer.propTypes = {
    data: PropTypes.oneOfType([PropTypes.object, PropTypes.array]).isRequired
};

export default DropdownContainer;
