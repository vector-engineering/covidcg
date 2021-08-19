import React from 'react';
import PropTypes from 'prop-types';
import { colors } from './VOCList';

import QuestionButton from '../Buttons/QuestionButton';
import {
  OrgLegendContainer,
  OrgLegendTitle,
  OrgItemContainer,
  OrgItem,
  OrgBadge,
  OrgName,
} from './OrgLegend.styles.js';

const Org = ({ name, color }) => {
  return (
    <OrgItem key={'org-' + name}>
      <OrgBadge color={color} />
      <OrgName>{name}</OrgName>
    </OrgItem>
  );
};
Org.propTypes = {
  name: PropTypes.string.isRequired,
  color: PropTypes.string.isRequired,
};

export const OrgLegend = () => {
  const orgs = [];
  Object.keys(colors).forEach((key) => {
    orgs.push(<Org key={'org-' + key} name={key} color={colors[key]} />);
  });

  return (
    <OrgLegendContainer>
      <OrgLegendTitle>Classified By</OrgLegendTitle>
      <OrgItemContainer>{orgs}</OrgItemContainer>
    </OrgLegendContainer>
  );
};
