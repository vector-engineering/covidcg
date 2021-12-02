import React from 'react';
import { observer } from 'mobx-react';
import useDimensions from 'react-use-dimensions';

// import ExternalLink from '../Common/ExternalLink';

import SurveillancePlot from '../Vega/SurveillancePlot';
import GlobalSeqPlot from '../Vega/GlobalSeqPlot';
import ExampleList from '../Example/ExampleList';

import {
  HomeTabContainer,
  HomeTabContent,
  // PubBanner,
  // CloseButton,
} from './HomeTab.styles';

const HomeTab = observer(() => {
  const [ref, { width }] = useDimensions();
  // const [showBanner, setShowBanner] = useState(true);

  return (
    <HomeTabContainer ref={ref}>
      {/* {showBanner && (
        <PubBanner>
          <p>
            COVID CG may not contain sequences submitted to GISAID after 2021-11-24. We are working on incorporating the missing data.
          </p>
          <CloseButton onClick={setShowBanner.bind(this, false)}>
            Dismiss
          </CloseButton>
        </PubBanner>
      )} */}
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
