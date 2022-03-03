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
            <p>Shown above are the dominant lineages per country within the last X months. Only the top 3 lineages globally are plotted above, other lineages are shown as grey. Countries without sufficient data are colored white.
</p>
            <MobileGlobalSeqPlot width={width}></MobileGlobalSeqPlot>
            <p>Shown above are the number of SARS-CoV-2 genomes (sequenced viruses—not cases!) since the start of the pandemic</p>
        </div>
    );
});

export default MobileHomePage;
