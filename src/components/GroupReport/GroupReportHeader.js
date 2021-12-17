import React from 'react';
import { observer } from 'mobx-react';

import VOCTable from './VOCList';

import { HeaderContainer, HeaderBoxRow } from './GroupReportHeader.styles';

const GroupReportHeader = observer(() => {
  return (
    <HeaderContainer>
      <HeaderBoxRow>
        <VOCTable />
      </HeaderBoxRow>
    </HeaderContainer>
  );
});

export default GroupReportHeader;
