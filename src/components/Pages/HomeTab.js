import React from 'react';
import { observer } from 'mobx-react';
import useDimensions from 'react-use-dimensions';

import SurveillancePlot from '../Vega/SurveillancePlot';
import GlobalSeqPlot from '../Vega/GlobalSeqPlot';
import ExampleList from '../Example/ExampleList';

import { HomeTabContainer, HomeTabContent } from './HomeTab.styles';

const HomeTab = observer(() => {
  const [ref, { width }] = useDimensions();
  return (
    <HomeTabContainer ref={ref}>
      <HomeTabContent>
        <SurveillancePlot width={width - 150} />
        <div style={{ height: '15px' }} />
        <GlobalSeqPlot width={width - 120} />
        <ExampleList />
      </HomeTabContent>
    </HomeTabContainer>
  );
});

export default HomeTab;
