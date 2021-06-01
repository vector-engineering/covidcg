import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';

import vocList from '../../../static_data/voc.json';

import {
  VOCListContainer,
  VOCListHeader,
  VOCListTitle,
  VOCItemListTitle,
  VOCItemList,
  VOCItemContainer,
  VOCItemName,
  VOCItemAlias,
} from './VOCList.styles';

const VOCItem = ({ name, alias }) => {
  return (
    <VOCItemContainer>
      <VOCItemName>{name}</VOCItemName>
      <VOCItemAlias>{alias}</VOCItemAlias>
    </VOCItemContainer>
  );
};
VOCItem.propTypes = {
  name: PropTypes.string,
  alias: PropTypes.string,
};

const VOCList = observer(() => {
  // Variants of Concern (VOC)
  const vocItems = [];
  vocList.list
    .filter((record) => {
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
  vocList.list
    .filter((record) => {
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
      <VOCItemListTitle>Variants of Concern</VOCItemListTitle>
      <VOCItemList>{vocItems}</VOCItemList>
      <VOCItemListTitle>Variants of Interest</VOCItemListTitle>
      <VOCItemList>{voiItems}</VOCItemList>
    </VOCListContainer>
  );
});

export default VOCList;
