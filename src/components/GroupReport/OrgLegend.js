import React from 'react';
import PropTypes from 'prop-types';
import ExternalLink from '../Common/ExternalLink';
import { colors } from './VOCList';

import {
  OrgLegendContainer,
  OrgLegendTitle,
  OrgItemContainer,
  OrgItem,
  OrgBadge,
  OrgName,
} from './OrgLegend.styles.js';

const urls = {
  WHO: 'https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/',
  CDC: 'https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-info.html',
  ECDC: 'https://www.ecdc.europa.eu/en/covid-19/variants-concern',
  PHE:
    'https://www.gov.uk/government/publications/covid-19-variants-genomically-confirmed-case-numbers',
};

const Org = ({ name, color }) => {
  return (
    <OrgItem key={'org-' + name}>
      <OrgBadge color={color} />
      <ExternalLink showIcon={false} href={urls[name]}>
        <OrgName>{name}</OrgName>
      </ExternalLink>
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
