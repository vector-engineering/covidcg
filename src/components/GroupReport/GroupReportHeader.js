import React from 'react';
import { observer } from 'mobx-react';

import GroupSearch from './GroupSearch';
import SelectedGroups from './SelectedGroups';
import VOCList from './VOCList';

import {
  HeaderContainer,
  Title,
  HeaderBoxRow,
} from './GroupReportHeader.styles';

const GroupReportHeader = observer(() => {
  return (
    <HeaderContainer>
      <Title>Lineage Report</Title>
      <HeaderBoxRow>
        <GroupSearch />
        <VOCList />
        <SelectedGroups />
      </HeaderBoxRow>
    </HeaderContainer>
  );
});

export default GroupReportHeader;
