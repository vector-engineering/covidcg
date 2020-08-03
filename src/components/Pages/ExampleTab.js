import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import styled from 'styled-components';

import { TABS } from '../../constants/UI';
import {
  GROUP_KEYS,
  DNA_OR_AA,
  COORDINATE_MODES,
} from '../../constants/config';

const ExampleTabContainer = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  padding: 20px;
`;

const ExampleHeader = styled.div`
  padding-left: 10px;
`;

const ExampleTitle = styled.h2``;

const ExampleList = styled.div`
  display: flex;
  flex-direction: row;
  flex-wrap: wrap;
`;
const ExampleItem = styled.a`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  width: 300px;
  height: 300px;
  margin: 10px;
  border: 1px solid #ccc;
  box-shadow: 0px 3px 4px #aaa;

  text-decoration: none;

  transition: 0.1s all ease-in-out;
  &:hover {
    border: 1px solid #00e;
  }
`;
const ExampleItemImage = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: center;

  img {
    width: auto;
    height: 240px;
  }
`;
const ExampleItemFooter = styled.div`
  height: 60px;
  padding: 10px;

  border-top: 1px solid #aaa;

  .example-item-title {
    font-size: 1.25em;
    font-weight: 700;
    color: #444;
  }
`;

const exampleItems = [
  {
    key: 'spike_d614g_na_eu',
    title: 'Rise of Spike D614G in Europe/North America',
    settings: {
      config: {
        groupKey: GROUP_KEYS.GROUP_SNV,
      },
      plotSettings: {},
      UI: {
        activeTab: TABS.TAB_GROUP,
      },
    },
  },
  // {
  //   key: 'b',
  //   title: 'Example B',
  // },
  // {
  //   key: 'c',
  //   title: 'Example C',
  // },
  // {
  //   key: 'd',
  //   title: 'Example D',
  // },
];

const ExampleTab = observer(() => {
  const onExampleClick = (key, e) => {
    e.preventDefault();
  };

  const exampleElements = [];
  exampleItems.forEach((exampleItem) => {
    exampleElements.push(
      <ExampleItem
        key={`example-item-${exampleItem.key}`}
        href="#"
        onClick={onExampleClick.bind(this, exampleItem.key)}
      >
        <ExampleItemImage>
          <img src="assets/images/cg_short_v13@4x_square.png" />
        </ExampleItemImage>
        <ExampleItemFooter>
          <span className="example-item-title">{exampleItem.title}</span>
        </ExampleItemFooter>
      </ExampleItem>
    );
  });

  return (
    <ExampleTabContainer>
      <ExampleHeader>
        <ExampleTitle>Example Analyses</ExampleTitle>
        <p>Get started with a bunch of example analyses</p>
      </ExampleHeader>
      <ExampleList>{exampleElements}</ExampleList>
    </ExampleTabContainer>
  );
});
ExampleTab.propTypes = {
  width: PropTypes.number,
};
ExampleTab.defaultProps = {
  width: 100,
};

export default ExampleTab;
