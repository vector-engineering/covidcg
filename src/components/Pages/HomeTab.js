import React, { useState } from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import useDimensions from 'react-use-dimensions';

import WalkthroughList from '../Example/WalkthroughList';
import SurveillancePlot from '../Viz/SurveillancePlot';
import GlobalSeqPlot from '../Viz/GlobalSeqPlot';
import ExampleList from '../Example/ExampleList';
import AcknowledgementFooter from '../Common/AcknowledgementFooter';

import { config } from '../../config';

import {
  HomeTabContainer,
  HomeTabContent,
  PubBanner,
  CloseButton,
} from './HomeTab.styles';

const HomeTab = observer(() => {
  const { UIStore } = useStores();
  const [ref, { width }] = useDimensions();
  const [showBanner, setShowBanner] = useState(config['show_home_banner']);

  return (
    <HomeTabContainer ref={ref}>
      {showBanner && UIStore.motd.length > 0 && (
        <PubBanner>
          <div dangerouslySetInnerHTML={{ __html: UIStore.motd }} />
          <CloseButton onClick={setShowBanner.bind(this, false)}>
            Dismiss
          </CloseButton>
        </PubBanner>
      )}
      <HomeTabContent>
        {config['show_walkthroughs'] && <WalkthroughList />}
        {config['show_surveillance'] && (
          <SurveillancePlot width={width - 150} />
        )}
        <div style={{ height: '15px' }} />
        {config['show_global_seq_plot'] && (
          <GlobalSeqPlot width={width - 120} />
        )}
        <ExampleList />
      </HomeTabContent>
      <AcknowledgementFooter />
    </HomeTabContainer>
  );
});

export default HomeTab;
