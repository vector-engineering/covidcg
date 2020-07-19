import React from 'react';
import Vega from './Vega';

export default function createClassFromSpec({ mode, spec }) {
  class FixedVegaChart extends React.PureComponent {
    static getSpec = function getSpec() {
      return spec;
    };

    render() {
      return <Vega mode={mode} spec={spec} {...this.props} />;
    }
  }

  return FixedVegaChart;
}
