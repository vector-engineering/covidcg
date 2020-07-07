import React from 'react';
import styled from 'styled-components';
import { Vega } from 'react-vega';
import { toJS } from 'mobx';
import AccordianWrapper from './AccordianWrapper';

const StyledTree = styled.div`
  .vega-bindings {
    display: flex;
    .vega-bind {
      margin-right: 12px;
    }
  }
`;

const VegaTree = ({ data, width }) => {
  const _data = [...toJS(data)];
  _data.push({
    id: 'root_node',
    name: 'root node',
    cases_percent: null,
    cases_sum: null,
    group: 'root_node',
  });

  console.log(width - 200);

  const vegaSpec = {
    $schema: 'https://vega.github.io/schema/vega/v5.json',
    description:
      'An example of Cartesian layouts for a node-link diagram of hierarchical data.',
    width: width - 200,
    height: 500,
    padding: 5,

    signals: [
      {
        name: 'labels',
        value: true,
        bind: { input: 'checkbox' },
      },
      {
        name: 'layout',
        value: 'tidy',
        bind: { input: 'radio', options: ['tidy', 'cluster'] },
      },
      {
        name: 'links',
        value: 'orthogonal',
        bind: {
          input: 'select',
          options: ['line', 'curve', 'diagonal', 'orthogonal'],
        },
      },
      {
        name: 'separation',
        value: true,
        bind: { input: 'checkbox' },
      },
    ],

    data: [
      {
        name: 'tree',
        values: _data,
        transform: [
          {
            type: 'stratify',
            key: 'id',
            parentKey: 'parent',
          },
          {
            type: 'tree',
            method: { signal: 'layout' },
            size: [{ signal: 'height' }, { signal: 'width' }],
            separation: { signal: 'separation' },
            as: ['y', 'x', 'depth', 'children'],
          },
        ],
      },
      {
        name: 'links',
        source: 'tree',
        transform: [
          { type: 'treelinks' },
          {
            type: 'linkpath',
            orient: 'horizontal',
            shape: { signal: 'links' },
          },
        ],
      },
    ],

    scales: [
      {
        name: 'color',
        type: 'linear',
        range: { scheme: 'magma' },
        domain: { data: 'tree', field: 'depth' },
        zero: true,
      },
    ],

    marks: [
      {
        type: 'path',
        from: { data: 'links' },
        encode: {
          update: {
            path: { field: 'path' },
            stroke: { value: '#ccc' },
          },
        },
      },
      {
        type: 'symbol',
        from: { data: 'tree' },
        encode: {
          enter: {
            size: { value: 100 },
            stroke: { value: '#fff' },
          },
          update: {
            x: { field: 'x' },
            y: { field: 'y' },
            fill: { scale: 'color', field: 'depth' },
          },
        },
      },
      {
        type: 'text',
        from: { data: 'tree' },
        encode: {
          enter: {
            text: { field: 'group' },
            fontSize: { value: 9 },
            baseline: { value: 'middle' },
          },
          update: {
            x: { field: 'x' },
            y: { field: 'y' },
            dx: { signal: 'datum.children ? -7 : 7' },
            align: { signal: "datum.children ? 'right' : 'left'" },
            opacity: { signal: 'labels ? 1 : 0' },
          },
        },
      },
    ],
  };

  console.log(_data);

  return (
    <AccordianWrapper
      defaultCollapsed={true}
      maxHeight={'1200px'}
      title={'tree'}
    >
      <StyledTree>
        <Vega spec={vegaSpec} signalListeners={{}} />
      </StyledTree>
    </AccordianWrapper>
  );
};

export default VegaTree;
