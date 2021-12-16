import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { toJS } from 'mobx';
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
import initialSpec from '../../vega_specs/num_seq_per_location_line.vg.json';

const PlotContainer = styled.div``;

const NumSeqPerLocationLine = observer(({ width }) => {
    const vegaRef = useRef();
    const {
        dataStore,
        configStore,
        UIStore,
    } = useStores();

    const handleHoverLocation = (...args) => {
        let hoverLocation;
        // Don't fire the action if there's no change
        if (args[1] == null) {
            hoverLocation = null;
        } else if (args[1]['location'] === undefined) {
            hoverLocation = null;
        } else {
            hoverLocation = args[1]['location'];
        }
        if (hoverLocation === configStore.hoverLocation) {
            return;
        }
        console.log('line hov', hoverLocation);
        configStore.updateHoverLocation(hoverLocation);
    };

    const processLocation = () => {
        let dat = toJS(dataStore.countsPerLocationDateMap);
        let dateTransformedData = [];
        // Get list of date for each key and iterate through all lists.
        // Potentially optimizable.
        dat.forEach((value, key) => {
            value.forEach((amount, date) => {
                // Convert date from unix time stamp to human readable, vega useable.
                dateTransformedData.push({ 'c': key, 'x': date, 'y': amount });
            });
        });
        dateTransformedData.sort(function (a, b) {
            let textA = a.c.toUpperCase();
            let textB = b.c.toUpperCase();
            return (textA > textB) ? -1 : (textA < textB) ? 1 : 0;
        });
        return dateTransformedData;
    };

    const [state, setState] = useState({
        data: {
            line_data: []
        },
        hoverLocation: null,
        spec: JSON.parse(JSON.stringify(initialSpec)),
        signalListeners: {
            hoverLocation: throttle(handleHoverLocation, 100)
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
                line_data: processLocation()
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
NumSeqPerLocationLine.propTypes = {
    width: PropTypes.number,
};
NumSeqPerLocationLine.defaultProps = {
    width: 100,
};

export default NumSeqPerLocationLine;
