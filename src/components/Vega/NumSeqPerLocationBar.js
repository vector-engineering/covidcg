import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { throttle, debounce } from '../../utils/func';

import EmptyPlot from '../Common/EmptyPlot';
import VegaEmbed from '../../react_vega/VegaEmbed';
import SkeletonElement from '../Common/SkeletonElement';
import { PlotOptions } from './Plot.styles';

import {
    ASYNC_STATES,
} from '../../constants/defs.json';
import initialSpec from '../../vega_specs/num_seq_per_location_bar.vg.json';

const PlotContainer = styled.div``;

const NumSeqPerLocationBar = observer(({ width }) => {
    const vegaRef = useRef();
    const {
        dataStore,
        configStore,
        UIStore,
    } = useStores();

    const handleHoverLocation = (...args) => {
        // Don't fire the action if there's no change
        let hoverLocation = args[1] === null ? null : args[1]['location'];
        if (hoverLocation === configStore.hoverLocation) {
            return;
        }
        configStore.updateHoverLocation(hoverLocation);
    };

    const processLocation = () => {
        let dat = dataStore.countsPerLocationMap;
        let locationData = [];
        for (var country in dat) {
            locationData.push({ 'c': country, 'amount': dat[country] });
        }
        return locationData;
    };

    const [state, setState] = useState({
        data: {
            bar_data: []
        },
        hoverLocation: null,
        spec: JSON.parse(JSON.stringify(initialSpec)),
        signalListeners: {
            hoverLocation: throttle(handleHoverLocation, 20)
        }
    });

    useEffect(() => {
        setState({
            ...state,
            hoverLocation: { location: configStore.hoverLocation },
        });
    }, [configStore.hoverLocation]);

    useEffect(() => {
        setState({
            ...state,
            data: {
                ...state.data
            },
        });
    }, [configStore.focusedLocations]);

    const refreshData = () => {
        if (UIStore.caseDataState !== ASYNC_STATES.SUCCEEDED) {
            return;
        }

        setState({
            ...state,
            data: {
                ...state.data,
                bar_data: processLocation()
            },
        });
    };

    // Refresh data on mount (i.e., tab change) or when data state changes
    useEffect(refreshData, [
        UIStore.caseDataState
    ]);
    useEffect(refreshData, []);

    if (UIStore.caseDataState === ASYNC_STATES.STARTED) {
        return (
            <div
                style={{
                    paddingTop: '12px',
                    paddingRight: '24px',
                    paddingLeft: '12px',
                    paddingBottom: '24px',
                }}
            >
                <SkeletonElement delay={2} height={100} />
            </div>
        );
    }

    if (configStore.selectedLocationNodes.length == 0) {
        return (
            <EmptyPlot height={100}>
                <p>
                    No locations selected. Please select one or more locations from the
                    sidebar, under &quot;Selected Locations&quot;, to compare counts of{' '}
                    <b>{configStore.getGroupLabel()}</b> between them.
                </p>
            </EmptyPlot>
        );
    }

    return (
        <PlotContainer>
            <PlotOptions>
                <div className="spacer" />
            </PlotOptions>
            <div style={{ width: `${width}` }}>
                <VegaEmbed
                    ref={vegaRef}
                    data={state.data}
                    spec={state.spec}
                    signals={{
                        hoverLocation: state.hoverLocation
                    }}
                    signalListeners={state.signalListeners}
                    dataListeners={state.dataListeners}
                    //width={width}
                    actions={false}
                />
            </div>
        </PlotContainer>
    );
});
NumSeqPerLocationBar.propTypes = {
    width: PropTypes.number,
};
NumSeqPerLocationBar.defaultProps = {
    width: 100,
};

export default NumSeqPerLocationBar;
