import React, { useState } from 'react';
import { observer } from 'mobx-react';
import useDimensions from 'react-use-dimensions';

import ExternalLink from '../Common/ExternalLink';

import MobileGlobalSeqPlot from '../Vega/MobileGlobalSeqPlot';
import CGLogo from '../../assets/images/cg_logo_v13.png';


import {
    MobileHomeContainer,
    MobileHomeContent,
    PubBanner,
    CloseButton,
    BannerLogo
} from './MobileHomePage.styles';

import SpinningGlobeSequences from '../Vega/SpinningGlobeSequences';

const MobileHomePage = observer(() => {
    const [ref, { width }] = useDimensions();
    const [showBanner, setShowBanner] = useState(true);


    return (
        <div>
            <BannerLogo src={CGLogo}></BannerLogo>
            <SpinningGlobeSequences width={width}></SpinningGlobeSequences>
            <MobileGlobalSeqPlot width={width}></MobileGlobalSeqPlot>
        </div>
    );
});

export default MobileHomePage;
