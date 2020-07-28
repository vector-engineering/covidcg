import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { toJS } from 'mobx';
import AccordionWrapper from '../Common/AccordionWrapper';
import { asyncStates } from '../../stores/uiStore';
import SkeletonElement from '../Common/SkeletonElement';
import { useStores } from '../../stores/connect';

import VegaEmbed from '../../react_vega/VegaEmbed';

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
  const { uiStore } = useStores();

  if (uiStore.caseDataState === asyncStates.STARTED) {
    return (
      <div
        style={{
          paddingTop: '0px',
          paddingRight: '24px',
          paddingLeft: '12px',
          paddingBottom: '0px',
        }}
      >
        <SkeletonElement delay={1} height={'50px'} />
      </div>
    );
  }

  const vegaSpec = {
    $schema: 'https://vega.github.io/schema/vega/v5.json',
    description:
      'An example of Cartesian layouts for a node-link diagram of hierarchical data.',
    width: width - 200,
    height: 500,
    padding: 5,

    signals: [
      {
        name: 'layout',
        value: 'tidy',
        bind: { input: 'radio', options: ['tidy', 'cluster'] },
      },
      {
        name: 'links',
        value: 'line',
        bind: {
          input: 'select',
          options: ['line', 'curve', 'diagonal', 'orthogonal'],
        },
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
            separation: true,
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
            strokeWidth: 12,
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
            fill: { field: 'color' },
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
            opacity: { signal: '1' },
          },
        },
      },
    ],
  };

  // console.log(_data);

  return (
    <AccordionWrapper
      defaultCollapsed={true}
      maxHeight={'1200px'}
      title={'tree'}
    >
      <StyledTree>
        <VegaEmbed spec={vegaSpec} signalListeners={{}} />
      </StyledTree>
    </AccordionWrapper>
  );
};

VegaTree.propTypes = {
  data: PropTypes.arrayOf(PropTypes.object).isRequired,
  width: PropTypes.number.isRequired,
};

export default VegaTree;
