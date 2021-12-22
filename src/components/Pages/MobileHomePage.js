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
        <MobileHomeContainer ref={ref}>
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
            <MobileHomeContent>
                <SpinningGlobeSequences></SpinningGlobeSequences>
            </MobileHomeContent>
        </MobileHomeContainer>
    );
});

export default MobileHomePage;
