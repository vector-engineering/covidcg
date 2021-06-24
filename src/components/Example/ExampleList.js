import React from 'react';
import { useStores } from '../../stores/connect';
import { observer } from 'mobx-react';

import ExternalLink from '../Common/ExternalLink';

import {
  ExampleListContainer,
  ExampleHeader,
  ExampleTitle,
  ExampleItem,
  ExampleItemList,
  ExampleItemImage,
  ExampleItemFooter,
} from './ExampleList.styles';

import TempImage from '../../assets/images/cg_short_v13@4x_square.png';

import { getExampleItems } from './examples';

const ExampleList = observer(() => {
  const { configStore, plotSettingsStore, UIStore, locationDataStore } =
    useStores();

  const selectTree = locationDataStore.selectTree;
  const exampleItems = getExampleItems({ selectTree });

  const onExampleClick = (title, e) => {
    e.preventDefault();

    const exampleItem = exampleItems.find((item) => item.title == title);
    // console.log(exampleItem);

    // Apply settings for each store
    Object.keys(exampleItem.settings).forEach((store) => {
      // Call the example action for each store
      if (store === 'config') {
        configStore.resetValues(exampleItem.settings[store]);
      } else if (store === 'plotSettings') {
        plotSettingsStore.resetValues(exampleItem.settings[store]);
      } else if (store === 'UI') {
        UIStore.resetValues(exampleItem.settings[store]);
      }
    });
  };

  const exampleElements = [];
  exampleItems.forEach((exampleItem) => {
    exampleElements.push(
      <ExampleItem
        key={`example-item-${exampleItem.title}`}
        href="#"
        onClick={onExampleClick.bind(this, exampleItem.title)}
      >
        <ExampleItemImage>
          <img
            src={
              exampleItem.image === undefined ? TempImage : exampleItem.image
            }
          />
        </ExampleItemImage>
        <ExampleItemFooter>
          <span className="example-item-title">{exampleItem.title}</span>
          <p className="example-item-description">{exampleItem.description}</p>
        </ExampleItemFooter>
      </ExampleItem>
    );
  });

  const renderExamples = () => {
    return exampleElements;
  };

  return (
    <ExampleListContainer>
      <ExampleHeader>
        <ExampleTitle>Example Analyses</ExampleTitle>
        <p>
          Use these example analyses to get started and explore the features of
          this application. If you would like to add an analysis to this list,
          please{' '}
          <ExternalLink href="https://github.com/vector-engineering/covidcg">
            submit a pull request on our GitHub
          </ExternalLink>
          , or contact us at{' '}
          <ExternalLink href="mailto:covidcg@broadinstitute.org">
            covidcg@broadinstitute.org
          </ExternalLink>
        </p>
      </ExampleHeader>
      <ExampleItemList>{renderExamples()}</ExampleItemList>
    </ExampleListContainer>
  );
});

export default ExampleList;
