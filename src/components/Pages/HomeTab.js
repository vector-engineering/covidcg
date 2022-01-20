import React, { useState } from 'react';
import { observer } from 'mobx-react';
import useDimensions from 'react-use-dimensions';

import ExternalLink from '../Common/ExternalLink';

import WalkthroughList from '../Example/WalkthroughList';
import SurveillancePlot from '../Vega/SurveillancePlot';
import GlobalSeqPlot from '../Vega/GlobalSeqPlot';
import ExampleList from '../Example/ExampleList';

import {
  HomeTabContainer,
  HomeTabContent,
  PubBanner,
  CloseButton,
} from './HomeTab.styles';

const HomeTab = observer(() => {
  const [ref, { width }] = useDimensions();
  const [showBanner, setShowBanner] = useState(true);

  return (
    <HomeTabContainer ref={ref}>
      {showBanner && (
        <PubBanner>
          <p>
            COVID CG is{' '}
            <ExternalLink href="https://doi.org/10.7554/eLife.63409">
              published in eLife
            </ExternalLink>
          </p>
          <CloseButton onClick={setShowBanner.bind(this, false)}>
            Dismiss
          </CloseButton>
        </PubBanner>
      )}
      <HomeTabContent>
        <WalkthroughList />
        <SurveillancePlot width={width - 150} />
        <div style={{ height: '15px' }} />
        <GlobalSeqPlot width={width - 120} />
        <ExampleList />
      </HomeTabContent>
    </HomeTabContainer>
  );
});

export default HomeTab;
