import React, { useState } from 'react';
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
  DropdownGrid,
  DropdownLabel,
} from './VOCList.styles';

const VOCItemDropdown = observer(({ data }) => {
  return (
    <DropdownGrid>
      <DropdownLabel>WHO Label: {data.who_label}</DropdownLabel>
      <DropdownLabel>Nextstrain: {data.nextstrain}</DropdownLabel>
      <DropdownLabel>Detected: {data.first_detection}</DropdownLabel>
    </DropdownGrid>
  );
});

const VOCItem = observer(({ data }) => {
  const { groupDataStore } = useStores();
  const [state, setState] = useState({
    open: false,
  });

  // See if checkbox should be checked
  let checked =
    groupDataStore.selectedGroups.indexOf(data.name) > -1 ? true : false;

  const checkBoxOnClick = (event) => {
    let selectedNodes = groupDataStore.selectedGroups;

    // If the checkbox is being checked, add to selectedNodes
    if (event.target.checked) {
      selectedNodes.push(data.name);
    } else {
      // If the check is being removed, remove from selectedNodes
      const index = selectedNodes.indexOf(data.name);
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

  const vocItemOnClick = () => {
    setState((state) => {
      return {
        open: !state.open,
      };
    });
  };

  return (
    <VOCItemContainer>
      {!state.open && (
        <VOCItemName onClick={vocItemOnClick}>▼ {data.name}</VOCItemName>
      )}
      {state.open && (
        <VOCItemName onClick={vocItemOnClick}>► {data.name}</VOCItemName>
      )}
      {state.open && <VOCItemDropdown data={data} />}
      {checked && (
        <input type="checkbox" onChange={checkBoxOnClick} data={data} checked />
      )}
      {!checked && (
        <input type="checkbox" onChange={checkBoxOnClick} data={data} />
      )}
    </VOCItemContainer>
  );
});
VOCItem.propTypes = {
  data: PropTypes.object,
};

const VOCList = observer(() => {
  // Variants of Concern (VOC)
  const vocItems = [];
  VOC_LIST.filter((record) => {
    return record.level === 'VOC';
  }).forEach((record) => {
    vocItems.push(<VOCItem key={`voc-item-${record.name}`} data={record} />);
  });

  // Variants of Interest (VOI)
  const voiItems = [];
  VOC_LIST.filter((record) => {
    return record.level === 'VOI';
  }).forEach((record) => {
    voiItems.push(<VOCItem key={`voi-item-${record.name}`} data={record} />);
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
