import React from 'react';
import PropTypes from 'prop-types';

import {
  LocationItemContainer,
  LocationItemLabel,
  DeselectButton,
} from './LocationItem.styles';

const LocationItem = ({ label, path, onDeselect }) => {
  const handleDeselect = (e) => {
    e.preventDefault();
    onDeselect(path);
  };

  return (
    <LocationItemContainer>
      <LocationItemLabel>{label}</LocationItemLabel>
      <DeselectButton title={'Deselect'} onClick={handleDeselect}>
        ×
      </DeselectButton>
    </LocationItemContainer>
  );
};
LocationItem.propTypes = {
  label: PropTypes.string.isRequired,
  path: PropTypes.string.isRequired,
  onDeselect: PropTypes.func.isRequired,
};

export default LocationItem;
