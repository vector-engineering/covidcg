import React from 'react';
import PropTypes from 'prop-types';

import { ListContainer, ItemContainer } from './StructureEntities.styles';

const StructureEntityItem = ({
  id,
  pdbx_description,
  type,
  poly_type,
  chains,
  checked,
  changeApplyHeatmap,
}) => {
  const chainsText = chains.join(', ');
  let typeText;
  if (poly_type !== null) {
    typeText = poly_type;
  } else {
    typeText = type;
  }

  return (
    <ItemContainer>
      <td>
        <input
          type="checkbox"
          checked={checked}
          onChange={changeApplyHeatmap.bind(this, id)}
        />
      </td>
      <td>{id}</td>
      <td>{pdbx_description}</td>
      <td>{typeText}</td>
      <td>{chainsText}</td>
    </ItemContainer>
  );
};
StructureEntityItem.propTypes = {
  id: PropTypes.string.isRequired,
  pdbx_description: PropTypes.string,
  type: PropTypes.string.isRequired,
  poly_type: PropTypes.string,
  seq: PropTypes.string,
  chains: PropTypes.arrayOf(PropTypes.string).isRequired,
  checked: PropTypes.bool,
  changeApplyHeatmap: PropTypes.func.isRequired,
};
StructureEntityItem.defaultProps = {
  poly_type: null,
  checked: true,
};

const StructureEntities = ({ entities, onChangeEntities }) => {
  const onChangeSingleEntityApplyHeatmap = (entityId) => {
    const newEntities = entities.map((e) => {
      e = Object.assign({}, e);
      if (e.id === entityId) {
        e.checked = !e.checked;
      }
      return e;
    });
    onChangeEntities(newEntities);
  };

  const entityItems = [];
  entities.forEach((entity) => {
    entityItems.push(
      <StructureEntityItem
        key={`structure-entity-item-${entity.id}`}
        changeApplyHeatmap={onChangeSingleEntityApplyHeatmap}
        {...entity}
      />
    );
  });

  return (
    <ListContainer>
      <table>
        <tbody>
          <tr>
            <th>Apply Heatmap</th>
            <th>ID</th>
            <th>Description</th>
            <th>Type</th>
            <th>Chains</th>
          </tr>
          {entityItems}
        </tbody>
      </table>
    </ListContainer>
  );
};
StructureEntities.propTypes = {
  entities: PropTypes.arrayOf(PropTypes.object).isRequired,
  onChangeEntities: PropTypes.func.isRequired,
};

export default StructureEntities;
