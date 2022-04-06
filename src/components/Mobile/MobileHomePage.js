import React, { useState } from 'react';
import { observer } from 'mobx-react';
import useDimensions from 'react-use-dimensions';

import ExternalLink from '../Common/ExternalLink';

import MobileGlobalSeqPlot from '../Vega/MobileGlobalSeqPlot';
import CGLogo from '../../assets/images/cg_logo_v13.png';


import {
    BannerLogo,
    MobileText
} from './MobileHomePage.styles';

import SpinningGlobeSequences from '../Vega/SpinningGlobeSequences';

const MobileHomePage = observer(() => {
    const [ref, { width }] = useDimensions();
    console.log(width);

    return (
        <div>
            <BannerLogo src={CGLogo}></BannerLogo>
            <SpinningGlobeSequences width={width}></SpinningGlobeSequences>
            <MobileText>Shown above are the dominant lineages per country within the last X months. Only the top 3 lineages globally are plotted above, other lineages are shown as grey. Countries without sufficient data are colored white.
</MobileText>
            <MobileGlobalSeqPlot width={width}></MobileGlobalSeqPlot>
            <MobileText>Shown above are the number of SARS-CoV-2 genomes (sequenced viruses—not cases!) since the start of the pandemic</MobileText>
        </div>
    );
});

export default MobileHomePage;
