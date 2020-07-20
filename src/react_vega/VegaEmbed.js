import React, {
  createRef,
  useLayoutEffect,
  useState,
  useRef,
  useEffect,
} from 'react';
import PropTypes from 'prop-types';
import vegaEmbed, { vega } from 'vega-embed';
import shallowEqual from './utils/shallowEqual';
import getUniqueFieldNames from './utils/getUniqueFieldNames';
import isFunction from './utils/isFunction';
import { NOOP } from './constants';
import computeSpecChanges from './utils/computeSpecChanges';
import equal from 'fast-deep-equal';

const VegaEmbed = ({
  className,
  spec,
  data,
  signals,
  signalListeners,
  dataListeners,
  style,
  width,
  height,
  onError,
  onNewView,
  ...options
}) => {
  const containerRef = createRef();
  let viewPromise;
  const [initialized, setInitialized] = useState(false);
  const [_view, setView] = useState(null);

  const prevOptionsRef = useRef(null);
  const prevSpecRef = useRef(null);
  const prevDataRef = useRef(null);
  const prevSignalListenersRef = useRef(null);
  const prevDataListenersRef = useRef(null);

  const handleError = (error) => {
    onError(error);
    // eslint-disable-next-line no-console
    console.warn('VEGA ERROR', error);
    return undefined;
  };

  const addSignalListenersToView = (view, signalListeners) => {
    if (!signalListeners) {
      return;
    }
    const signalNames = Object.keys(signalListeners);
    signalNames.forEach((signalName) => {
      try {
        view.addSignalListener(signalName, signalListeners[signalName]);
      } catch (error) {
        // eslint-disable-next-line no-console
        console.warn('Cannot add invalid signal listener.', error);
      }
    });
    return signalNames.length > 0;
  };

  const removeSignalListenersFromView = (view, signalListeners) => {
    if (!signalListeners) {
      return;
    }
    const signalNames = Object.keys(signalListeners);
    signalNames.forEach((signalName) => {
      try {
        view.removeSignalListener(signalName, signalListeners[signalName]);
      } catch (error) {
        // eslint-disable-next-line no-console
        console.warn('Cannot remove invalid signal listener.', error);
      }
    });
    return signalNames.length > 0;
  };

  const addDataListenersToView = (view, dataListeners) => {
    if (!dataListeners) {
      return;
    }
    const dataNames = Object.keys(dataListeners);
    dataNames.forEach((dataName) => {
      try {
        view.addDataListener(dataName, dataListeners[dataName]);
      } catch (error) {
        console.warn('Cannot add invalid data listener.', error);
      }
    });
    return dataNames.length > 0;
  };

  const removeDataListenersFromView = (view, dataListeners) => {
    if (!dataListeners) {
      return;
    }
    const dataNames = Object.keys(dataListeners);
    dataNames.forEach((dataName) => {
      try {
        view.removeDataListener(dataName, dataListeners[dataName]);
      } catch (error) {
        console.warn('Cannot remove invalid data listener.', error);
      }
    });
    return dataNames.length > 0;
  };

  const combineSpecWithDimension = ({ spec, width, height }) => {
    if (typeof width !== 'undefined' && typeof height !== 'undefined') {
      return { ...spec, width, height };
    }
    if (typeof width !== 'undefined') {
      return { ...spec, width };
    }
    if (typeof height !== 'undefined') {
      return { ...spec, height };
    }
    return spec;
  };

  const createView = () => {
    // console.log('CREATE VIEW');
    if (containerRef.current) {
      const finalSpec = combineSpecWithDimension({ spec, width, height });
      viewPromise = vegaEmbed(containerRef.current, finalSpec, options)
        .then(({ view }) => {
          let hasListeners =
            addSignalListenersToView(view, signalListeners) +
            addDataListenersToView(view, dataListeners);
          if (hasListeners) {
            view.run();
          }
          setView(view);
          return { view };
        })
        .catch(handleError);

      if (onNewView) {
        modifyView(onNewView);
      }
    }
  };

  const modifyView = (action) => {
    // console.log('MODIFY VIEW', viewPromise);
    if (viewPromise) {
      viewPromise
        .then(({ view }) => {
          if (view) {
            action(view);
          }
          return { view };
        })
        .catch(handleError);
    }
  };

  const clearView = () => {
    // console.log('CLEAR VIEW');
    modifyView((view) => {
      view.finalize();
    });
    viewPromise = undefined;

    return this;
  };

  useEffect(() => {
    prevSpecRef.current = spec;
    prevOptionsRef.current = options;
    prevSignalListenersRef.current = signalListeners;
    prevDataListenersRef.current = dataListeners;
  });

  // Listen to dataset changes and update data
  useEffect(() => {
    // console.log('NEW DATA');
    // console.log(data);
    // console.log(viewPromise);

    // I don't know what I'm doing. but the promise is defined
    // during initialization, so this runs on hard re-render
    if (viewPromise && data && Object.keys(data).length > 0) {
      modifyView((view) => {
        Object.keys(data).forEach((name) => {
          if (data[name]) {
            if (isFunction(data[name])) {
              data[name](view.data(name));
            } else {
              view.change(
                name,
                vega
                  .changeset()
                  .remove(() => true)
                  // Do another deep copy of the data... since Vega will insert
                  // additional keys during its processing and we want an unmodified
                  // version to do diffs with
                  .insert(JSON.parse(JSON.stringify(data[name])))
              );
            }
          }
        });
        view.resize().run();
      });
      // And this will run on new props (soft re-render)
    } else if (_view && data) {
      let prevData = prevDataRef.current;

      Object.keys(data).forEach((name) => {
        // console.log(name, data[name], prevData[name]);

        if (
          !Object.prototype.hasOwnProperty.call(prevData, name) ||
          !equal(data[name], prevData[name])
        ) {
          console.log('Changing datasets', name);
          if (isFunction(data[name])) {
            data[name](_view.data(name));
          } else {
            _view.change(
              name,
              vega
                .changeset()
                .remove(() => true)
                .insert(JSON.parse(JSON.stringify(data[name])))
            );
          }
        }
      });
    }

    prevDataRef.current = data;
  }, [data]);

  // Update Vega dimensions without redrawing the whole thing
  useEffect(() => {
    // console.log('change width/height');
    if (_view) {
      _view.width(width);
      _view.height(height);
    }
  }, [width, height]);

  // Listen to changes in signals passed via. props
  useEffect(() => {
    // console.log('Passed signals:', signals);
    if (_view) {
      Object.keys(signals).forEach((signalName) => {
        let currentSignalVal = _view.signal(signalName);
        // console.log('current', currentSignalVal, 'new', signals[signalName]);
        // Only update if the signal is different
        if (currentSignalVal != signals[signalName]) {
          _view.signal(signalName, { group: signals[signalName] });
        }
      });
    }
  }, [signals]);

  // Listen to changes in signalListeners
  const updateSignalListeners = () => {
    const newSignalListeners = signalListeners;
    const oldSignalListeners = prevSignalListenersRef.current;

    const areSignalListenersChanged = !shallowEqual(
      newSignalListeners,
      oldSignalListeners
    );

    if (_view && areSignalListenersChanged) {
      removeSignalListenersFromView(_view, oldSignalListeners);
      addSignalListenersToView(_view, newSignalListeners);
      _view.runAsync();
    }
  };
  useEffect(() => {
    updateSignalListeners();
  }, [signalListeners]);

  // Listen to changes in dataListeners
  const updateDataListeners = () => {
    const newDataListeners = dataListeners;
    const oldDataListeners = prevDataListenersRef.current;

    const areDataListenersChanged = !shallowEqual(
      newDataListeners,
      oldDataListeners
    );

    if (_view && areDataListenersChanged) {
      removeDataListenersFromView(_view, oldDataListeners);
      addDataListenersToView(_view, newDataListeners);
      _view.runAsync();
    }
  };
  useEffect(() => {
    updateDataListeners();
  }, [dataListeners]);

  const recreateView = () => {
    clearView();
    createView();
    updateSignalListeners();
    updateDataListeners();
  };

  useLayoutEffect(() => {
    if (!initialized) {
      // console.log('initializing');
      createView();
      setInitialized(true);
      return;
    }

    const prevOptions = prevOptionsRef.current;
    const fieldSet = getUniqueFieldNames([options, prevOptions]);
    // console.log(options, prevOptions, fieldSet);
    const prevSpec = prevSpecRef.current;

    // Only create a new view if necessary
    if (Array.from(fieldSet).some((f) => options[f] !== prevOptions[f])) {
      recreateView();
    } else {
      const specChanges = computeSpecChanges(spec, prevSpec);
      // console.log(spec, prevSpec);
      // console.log('spec changes', specChanges);

      if (specChanges && specChanges.isExpensive) {
        recreateView();
      }
    }

    // Specify how to clean up after this effect:
    // This function is called when the component unmounts
    return function cleanup() {
      clearView();
    };
  }, [spec, options]);

  return <div ref={containerRef} className={className} style={style} />;
};

VegaEmbed.propTypes = {
  mode: PropTypes.oneOf(['vega', 'vega-lite']),
  height: PropTypes.number,
  width: PropTypes.number,

  className: PropTypes.string,
  spec: PropTypes.object.isRequired,
  data: PropTypes.object.isRequired,
  signals: PropTypes.object,
  signalListeners: PropTypes.object, // key -> value (function)
  dataListeners: PropTypes.object, // key -> value (function)
  style: PropTypes.object,
  onNewView: PropTypes.func,
  onError: PropTypes.func,
};
VegaEmbed.defaultProps = {
  mode: 'vega',
  className: 'vega-embed',
  signals: {},
  signalListeners: {},
  dataListeners: {},
  style: {},
  onNewView: NOOP,
  onError: (error) => {
    console.error(error);
    return null;
  },
};

export default VegaEmbed;
