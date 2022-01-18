import React, { useState } from 'react';
import { observer } from 'mobx-react';
import useDimensions from 'react-use-dimensions';

import ExternalLink from '../Common/ExternalLink';

import MobileGlobalSeqPlot from '../Vega/MobileGlobalSeqPlot';

import {
    MobileHomeContainer,
    MobileHomeContent,
    PubBanner,
    CloseButton,
} from './MobileHomePage.styles';

import SpinningGlobeSequences from '../Vega/SpinningGlobeSequences';

const MobileHomePage = observer(() => {
    const [ref, { width }] = useDimensions();
    const [showBanner, setShowBanner] = useState(true);


    return (
        <MobileHomeContainer>
            <SpinningGlobeSequences width={400}></SpinningGlobeSequences>
            <MobileGlobalSeqPlot width={400}></MobileGlobalSeqPlot>
        </MobileHomeContainer>
    );
});

export default MobileHomePage;
