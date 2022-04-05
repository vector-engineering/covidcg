import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { asyncDataStoreInstance } from '../App';

import {
  VOCTableContainer,
  VOCTableRow,
  VOCTableHeader,
  VOCGridTitle,
  VOCTableContent,
  VOCItemContainer,
  VOCItemName,
} from './VOCList.styles';

const RSVGenotypeItem = observer(({ genotype }) => {
  const { groupDataStore } = useStores();
  // See if container should be selected
  let selected =
    groupDataStore.selectedGroups.indexOf(genotype) > -1 ? true : false;

  const onClick = (event) => {
    let selectedNodes = groupDataStore.selectedGroups;

    // If the container is selected, add to selectedNodes
    if (!event.target.selected) {
      selectedNodes.push(genotype);
    } else {
      // If the container is being deselected, remove from selectedNodes
      const index = selectedNodes.indexOf(genotype);
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
      <VOCItemName selected={selected}>{genotype}</VOCItemName>
    </VOCItemContainer>
  );
});
RSVGenotypeItem.propTypes = {
  genotype: PropTypes.string.isRequired,
};

const RSVGenotypeList = observer(() => {
  const { configStore } = useStores();
  const genotypes =
    asyncDataStoreInstance.data.genotypesBySubtype[
      configStore.selectedReference
    ];

  const items = [];

  genotypes.forEach((genotype) => {
    items.push(<RSVGenotypeItem genotype={genotype} key={genotype} />);
  });

  return (
    <VOCTableContainer>
      <VOCTableRow>
        <div>
          <VOCTableHeader>
            <VOCGridTitle>
              Genotypes for RSV-{configStore.selectedReference}
            </VOCGridTitle>
          </VOCTableHeader>
          <VOCTableContent>{items}</VOCTableContent>
        </div>
      </VOCTableRow>
    </VOCTableContainer>
  );
});

export default RSVGenotypeList;
