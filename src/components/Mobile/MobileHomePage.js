import React, { useState } from 'react';
import { observer } from 'mobx-react';
import useDimensions from 'react-use-dimensions';

import ExternalLink from '../Common/ExternalLink';

import SurveillancePlot from '../Vega/SurveillancePlot';
import GlobalSeqPlot from '../Vega/GlobalSeqPlot';
import ExampleList from '../Example/ExampleList';

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
        <SpinningGlobeSequences width={500}></SpinningGlobeSequences>
    );
});

export default MobileHomePage;
