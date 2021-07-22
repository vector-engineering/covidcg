import React, { useState } from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import VOC_LIST from '../../../static_data/vocs.json';

import {
  VOCListContainer,
  VOCGridTitle,
  GridItem,
  VOCItemContainer,
  VOCItemGrid,
  VOCItemName,
  VOCBadgeContainer,
  VOCBadge,
} from './VOCList.styles';

const colors = { WHO: '#88CCEE', CDC: '#DDCC77', ECDC: '#AA4499' };
const coords = { WHO: [1, 1], CDC: [1, 2], ECDC: [2, 1], EMPTY: [2, 2] };

const VOCItem = observer(({ name, orgArr }) => {
  const { groupDataStore } = useStores();
  const badges = [];

  // To display borders, four badges must be displayed
  // Loop through orgArr and if real, display colored badge, else display white
  Object.keys(coords).forEach((item) => {
    badges.push(
      <VOCBadge
        color={orgArr.includes(item) ? colors[item] : 'white'}
        row={coords[item][0]}
        col={coords[item][1]}
      />
    );
  });

  // See if container should be selected
  let selected =
    groupDataStore.selectedGroups.indexOf(name) > -1 ? true : false;

  const onClick = (event) => {
    let selectedNodes = groupDataStore.selectedGroups;

    // If the container is selected, add to selectedNodes
    if (!event.target.selected) {
      selectedNodes.push(name);
    } else {
      // If the container is being deselected, remove from selectedNodes
      const index = selectedNodes.indexOf(name);
      if (index > -1) {
        selectedNodes.splice(index, 1);
      }
    }

    // Update selectedNodes
    groupDataStore.updateSelectedGroups(
      selectedNodes.map((node) => {
        return node;
      })
    );
  };

  return (
    <VOCItemContainer onClick={onClick} selected={selected}>
      <VOCBadgeContainer>{badges}</VOCBadgeContainer>
      <VOCItemName selected={selected}>{name}</VOCItemName>
    </VOCItemContainer>
  );
});
VOCItem.propTypes = {
  name: PropTypes.string.isRequired,
  orgArr: PropTypes.arrayOf(PropTypes.string).isRequired,
};

const VOCList = observer(() => {
  const filterItems = (level) => {
    // VOC_LIST is an object with lineages as its keys
    // Each lineage is an object with organizations (WHO, CDC, etc) as its keys
    // Each organization's value is the level (VOC, VOI, Other) assigned to the lineage

    // {lineage: {
    //            org1: level,
    //            org2: level, ...
    //           }
    // }

    const items = [];
    Object.keys(VOC_LIST)
      .filter((lineage_name) => {
        return Object.values(VOC_LIST[lineage_name]).some((v) => v === level);
      })
      .forEach((lineage_name) => {
        let org_dict = VOC_LIST[lineage_name];
        let orgArr = Object.keys(org_dict).filter((org) => {
          return org_dict[org] === level;
        });

        items.push(
          <VOCItem
            key={`${level}-item-${lineage_name}`}
            name={lineage_name}
            orgArr={orgArr}
          />
        );
      });
    return items;
  };

  // Variants of Concern (VOC)
  const vocItems = filterItems('VOC');

  // Variants of Interest (VOI)
  const voiItems = filterItems('VOI');

  return (
    <VOCListContainer>
      <VOCItemGrid>
        <VOCGridTitle>Variants of Concern</VOCGridTitle>
        <VOCGridTitle>Variants of Interest</VOCGridTitle>
        <GridItem>{vocItems}</GridItem>
        <GridItem>{voiItems}</GridItem>
      </VOCItemGrid>
    </VOCListContainer>
  );
});

export default VOCList;