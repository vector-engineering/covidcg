import React from 'react';
import { colors } from './VOCList';

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
