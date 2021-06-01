import React from 'react';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';
import useDimensions from 'react-use-dimensions';

import GroupTreePlot from '../Vega/GroupTreePlot';
import GroupReportHeader from '../GroupReport/GroupReportHeader';

const GroupReportTabContainer = styled.div`
  display: flex;
  flex-direction: row;
`;
const GroupTreePlotContainer = styled.div``;

const MainContainer = styled.div``;

const GroupReportTab = observer(() => {
  const { configStore } = useStores();
  const [ref, { width }] = useDimensions();

  return (
    <GroupReportTabContainer ref={ref}>
      <GroupTreePlotContainer>
        <GroupTreePlot width={300} />
      </GroupTreePlotContainer>
      <MainContainer>
        <GroupReportHeader />
      </MainContainer>
    </GroupReportTabContainer>
  );
});

export default GroupReportTab;
