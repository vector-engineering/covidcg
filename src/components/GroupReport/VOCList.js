import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import VOC_LIST from '../../../static_data/vocs.json';

import {
  VOCListContainer,
  VOCListHeader,
  VOCListTitle,
  VOCTitle,
  VOITitle,
  VOCGridItem,
  VOIGridItem,
  VOCItemContainer,
  VOCItemGrid,
  VOCItemName,
  VOCItemAlias,
} from './VOCList.styles';

const VOCItem = observer(({ name, alias }) => {
  const { groupDataStore } = useStores();

  const vocItemOnClick = () => {
    let selectedNodes = groupDataStore.selectedGroups;

    let inSelectedNodes = false;

    // Check if the clicked node is selected
    for (let i = 0; i < selectedNodes.length; i++) {
      const node = selectedNodes[i];
      if (node.localeCompare(name) === 0) {
        // If the clicked variant is already selected, unselect it
        selectedNodes.splice(i, 1);
        inSelectedNodes = true;
        break;
      }
    }

    // If the clicked variant is not selected, select it
    if (!inSelectedNodes) {
      selectedNodes.push(name);
    }

    groupDataStore.updateSelectedGroups(
      selectedNodes.map((node) => {
        return node;
      })
    );
  };

  if (alias) {
    return (
      <VOCItemContainer onClick={vocItemOnClick}>
        <VOCItemName>{name} - </VOCItemName>
        <VOCItemAlias>{alias}</VOCItemAlias>
      </VOCItemContainer>
    );
  } else {
    return (
      <VOCItemContainer onClick={vocItemOnClick}>
        <VOCItemName>{name}</VOCItemName>
      </VOCItemContainer>
    );
  }
});
VOCItem.propTypes = {
  name: PropTypes.string,
  alias: PropTypes.string,
};

const VOCList = observer(() => {
  // Variants of Concern (VOC)
  const vocItems = [];
  VOC_LIST.filter((record) => {
      return record.level === 'VOC';
    })
    .forEach((record) => {
      vocItems.push(
        <VOCItem
          key={`voc-item-${record.name}`}
          name={record.name}
          alias={record.alias}
        />
      );
    });

  // Variants of Interest (VOI)
  const voiItems = [];
  VOC_LIST.filter((record) => {
      return record.level === 'VOI';
    })
    .forEach((record) => {
      voiItems.push(
        <VOCItem
          key={`voi-item-${record.name}`}
          name={record.name}
          alias={record.alias}
        />
      );
    });

  return (
    <VOCListContainer>
      <VOCListHeader>
        <VOCListTitle>Select Notable Variants</VOCListTitle>
      </VOCListHeader>
      <VOCItemGrid>
        <VOCTitle>Variants of Concern</VOCTitle>
        <VOITitle>Variants of Interest</VOITitle>
        <VOCGridItem>{vocItems}</VOCGridItem>
        <VOIGridItem>{voiItems}</VOIGridItem>
      </VOCItemGrid>
    </VOCListContainer>
  );
});

export default VOCList;
