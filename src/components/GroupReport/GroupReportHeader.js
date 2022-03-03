import React from 'react';
import { observer } from 'mobx-react';

import { config } from '../../config';

import VOCTable from './VOCList';
import RSVGenotypeList from './RSVGenotypeList';

import { HeaderContainer, HeaderBoxRow } from './GroupReportHeader.styles';

const GroupReportHeader = observer(() => {
  return (
    <HeaderContainer>
      <HeaderBoxRow>
        {config.virus === 'sars2' && <VOCTable />}
        {config.virus === 'rsv' && <RSVGenotypeList />}
      </HeaderBoxRow>
    </HeaderContainer>
  );
});

export default GroupReportHeader;
