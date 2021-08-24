import React from 'react';
import { observer } from 'mobx-react';

import VOCList from './VOCList';

import { HeaderContainer, HeaderBoxRow } from './GroupReportHeader.styles';

const GroupReportHeader = observer(() => {
  return (
    <HeaderContainer>
      <HeaderBoxRow>
        <VOCList />
      </HeaderBoxRow>
    </HeaderContainer>
  );
});

export default GroupReportHeader;
