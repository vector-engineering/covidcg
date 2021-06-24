import React from 'react';
import { observer } from 'mobx-react';

import GroupSearch from './GroupSearch';
import VOCList from './VOCList';

import { HeaderContainer, HeaderBoxRow } from './GroupReportHeader.styles';

const GroupReportHeader = observer(() => {
  return (
    <HeaderContainer>
      <HeaderBoxRow>
        <GroupSearch />
        <VOCList />
      </HeaderBoxRow>
    </HeaderContainer>
  );
});

export default GroupReportHeader;
