import React, { useState } from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { asyncDataStoreInstance } from '../App';

import { OrgLegend } from './OrgLegend';

import {
  VOCTableContainer,
  VOCTableRow,
  VOCTableHeader,
  VOCTableToggle,
  VOCTableContent,
  VOCGridTitle,
  VOCItemContainer,
  VOCBadgeContainer,
  VOCBadge,
} from './VOCList.styles';

import { LineageName } from '../Common/LineageName';

export const colors = {
  WHO: '#88CCEE',
  CDC: '#DDCC77',
  ECDC: '#AA4499',
  PHE: '#004488',
};
const coords = { WHO: [1, 1], CDC: [1, 2], ECDC: [2, 1], PHE: [2, 2] };

const VOCItem = observer(({ name, orgArr }) => {
  const { groupDataStore } = useStores();
  const badges = [];

  // To display borders, four badges must be displayed
  // Loop through orgArr and if real, display colored badge, else display white
  Object.keys(coords).forEach((item, index) => {
    badges.push(
      <VOCBadge
        color={orgArr.includes(item) ? colors[item] : 'white'}
        row={coords[item][0]}
        col={coords[item][1]}
        key={`${item}-badge-${index}`}
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
      <LineageName name={name} selected={selected} />
    </VOCItemContainer>
  );
});
VOCItem.propTypes = {
  name: PropTypes.string.isRequired,
  orgArr: PropTypes.arrayOf(PropTypes.string).isRequired,
};

const VOCTable = observer(() => {
  const filterItems = (level) => {
    // VOC_LIST is an object with lineages as its keys
    // Each lineage is an object with organizations (WHO, CDC, etc) as its keys
    // Each organization's value is the level (VOC, VOI, Other) assigned to the lineage

    // {lineage: {
    //            org1: level,
    //            org2: level, ...
    //           }
    // }

    const VOC_LIST = asyncDataStoreInstance.data.vocs;

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
  // The PHE calls this same level (one below VOC) VUI
  const voiItems = filterItems('VOI');

  // Other Variants
  const otherItems = filterItems('Other');

  const [state, setState] = useState({
    vocExpanded: true,
    voiExpanded: true,
    otherExpanded: false,
  });

  const toggleExpand = (event) => {
    const category = event.target.getAttribute('name');

    if (category === 'VOC') {
      setState({ ...state, vocExpanded: !state.vocExpanded });
    } else if (category === 'VOI') {
      setState({ ...state, voiExpanded: !state.voiExpanded });
    } else if (category === 'Other') {
      setState({ ...state, otherExpanded: !state.otherExpanded });
    }
  };

  return (
    <VOCTableContainer>
      <OrgLegend />
      <VOCTableRow>
        <div>
          <VOCTableHeader>
            <VOCGridTitle>Variants of Concern *</VOCGridTitle>
            <VOCTableToggle
              onClick={toggleExpand}
              name="VOC"
              expanded={state.vocExpanded}
            />
          </VOCTableHeader>
          {state.vocExpanded && <VOCTableContent>{vocItems}</VOCTableContent>}
          {!state.vocExpanded && <p>{vocItems.length} VOCs hidden...</p>}
        </div>
      </VOCTableRow>
      <VOCTableRow>
        <div>
          <VOCTableHeader>
            <VOCGridTitle>Variants of Interest</VOCGridTitle>
            <VOCTableToggle
              onClick={toggleExpand}
              name="VOI"
              expanded={state.voiExpanded}
            />
          </VOCTableHeader>
          {state.voiExpanded && <VOCTableContent>{voiItems}</VOCTableContent>}
          {!state.voiExpanded && <p>{voiItems.length} VOIs hidden...</p>}
        </div>
      </VOCTableRow>
      <VOCTableRow>
        <div>
          <VOCTableHeader>
            <VOCGridTitle>Other Variants Being Monitored</VOCGridTitle>
            <VOCTableToggle
              onClick={toggleExpand}
              name="Other"
              expanded={state.otherExpanded}
            />
          </VOCTableHeader>
          {state.otherExpanded && (
            <VOCTableContent>{otherItems}</VOCTableContent>
          )}
          {!state.otherExpanded && (
            <p>{otherItems.length} variants hidden...</p>
          )}
        </div>
      </VOCTableRow>
      <p>* The WHO, CDC, and ECDC classify all descendents of VOCs as VOCs</p>
    </VOCTableContainer>
  );
});

export default VOCTable;
