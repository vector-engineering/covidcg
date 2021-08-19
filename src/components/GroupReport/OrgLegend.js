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
      <QuestionButton
        data-tip={`
          <div>
          <p>Displayed below are the PANGO lineages of variants being specifically monitored by at least one of these organizations.</p>
          <p>Variants of Concern are being most closely monitored, followed by Variants of Interest (called Variants Under Investigation by the UK PHE).</p>
          <p>Any other, organization-specific classifications are collectively grouped under 'Other Variants Being Monitored'.</p>
          <p>If a square next to a lineage is colored in, that organization has given the lineage whichever classification it is below.</p>
          <p>(i.e. B.1.617.2 has been labeled a Variant of Concern by the WHO, CDC, ECDC, and PHE).</p>
          </div>
        `}
        data-html="true"
        data-for="main-tooltip"
      />
    </OrgLegendContainer>
  );
};
