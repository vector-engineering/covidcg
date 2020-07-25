import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { observer } from 'mobx-react';

const LocationTabContainer = styled.div``;

const LocationTab = observer(({ width }) => {
  return <LocationTabContainer></LocationTabContainer>;
});
LocationTab.propTypes = {
  width: PropTypes.number,
};
LocationTab.defaultProps = {
  width: 100,
};

export default LocationTab;
