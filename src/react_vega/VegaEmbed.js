import React, {
  useLayoutEffect,
  useState,
  useRef,
  useEffect,
  forwardRef,
  useImperativeHandle,
} from 'react';
import PropTypes from 'prop-types';
import vegaEmbed, { vega } from 'vega-embed';
import shallowEqual from './utils/shallowEqual';
import getUniqueFieldNames from './utils/getUniqueFieldNames';
import isFunction from './utils/isFunction';
import computeSpecChanges from './utils/computeSpecChanges';
import equal from 'fast-deep-equal';

const VegaEmbed = forwardRef(
  (
    {
      className,
      spec,
      data,
      recreateOnDatasets,
      signals,
      cheapSignals,
      updateDataSignals,
      signalListeners,
      dataListeners,
      style,
      width,
      height,
      onError,
      onComplete,
      onNewView,
      ...options
    },
    ref
  ) => {
    const containerRef = useRef();

    const [state, setState] = useState({
      viewPromise: undefined,
      spec: spec,
    });

    const prevOptionsRef = useRef(options);
    const prevDataRef = useRef(data);
    const prevSignalListenersRef = useRef(signalListeners);
    const prevDataListenersRef = useRef(dataListeners);

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

    const createView = () => {
      if (containerRef.current) {
        // console.log('CREATE VIEW');
        // vega-embed options: https://github.com/vega/vega-embed#options
        let viewPromise = vegaEmbed(containerRef.current, spec, options)
          .then(({ view }) => {
            // Add signals
            Object.keys(signals).forEach((signalName) => {
              view.signal(signalName, signals[signalName]);
            });

            return view;
          })
          .then((view) => {
            // Unset previous data to force push the data into the new view
            prevDataRef.current = {};
            updateData(view);

            return view;
          })
          .then((view) => {
            addSignalListenersToView(view, signalListeners);
            addDataListenersToView(view, dataListeners);

            if (onComplete) {
              view.runAfter(onComplete);
            }

            if (width !== undefined) {
              view.width(width);
            }
            if (height !== undefined) {
              view.height(height);
            }
            view.resize().run();
            return view;
          })
          .catch(handleError);

        setState({ ...state, viewPromise, spec: spec });

        if (onNewView !== null) {
          modifyView(onNewView);
        }
      }
    };

    const modifyView = (action) => {
      // console.log('MODIFY VIEW');
      if (state.viewPromise !== undefined) {
        state.viewPromise
          .then((view) => {
            if (view) {
              action(view);
            }
            return view;
          })
          .catch(handleError);
      }
    };

    const clearView = () => {
      modifyView((view) => {
        // console.log('CLEAR VIEW');
        view.finalize();
        // setState({ ...state, viewPromise: undefined });

        prevDataRef.current = {};
        prevSignalListenersRef.current = {};
        prevDataListenersRef.current = {};
      });
    };

    const updateData = (view) => {
      if (view && data) {
        let prevData = prevDataRef.current;
        let changed = false;

        Object.keys(data).forEach((name) => {
          // console.log(name, data[name], prevData[name]);
          // if (data[name] !== prevData[name]) {
          //   console.log('datasets changed:', name);
          //   console.log(data[name]);
          // }

          if (
            !Object.prototype.hasOwnProperty.call(prevData, name) ||
            data[name] !== prevData[name]
          ) {
            console.debug('VEGA Changing datasets', name);
            changed = true;
            if (isFunction(data[name])) {
              data[name](view.data(name));
            } else {
              view.change(
                name,
                vega
                  .changeset()
                  .remove(() => true)
                  .insert(JSON.parse(JSON.stringify(data[name])))
              );
            }
          }
        });

        if (changed) {
          view.resize().run();
        }
      }
      prevDataRef.current = data;
    };

    // Listen to changes in signalListeners
    const updateSignalListeners = () => {
      const newSignalListeners = signalListeners;
      const oldSignalListeners = prevSignalListenersRef.current;

      const areSignalListenersChanged = !shallowEqual(
        newSignalListeners,
        oldSignalListeners
      );

      if (areSignalListenersChanged) {
        modifyView((view) => {
          removeSignalListenersFromView(view, oldSignalListeners);
          addSignalListenersToView(view, newSignalListeners);
        });
      }
      prevSignalListenersRef.current = signalListeners;
    };

    // Listen to changes in dataListeners
    const updateDataListeners = () => {
      const newDataListeners = dataListeners;
      const oldDataListeners = prevDataListenersRef.current;

      const areDataListenersChanged = !shallowEqual(
        newDataListeners,
        oldDataListeners
      );

      if (areDataListenersChanged) {
        modifyView((view) => {
          removeDataListenersFromView(view, oldDataListeners);
          addDataListenersToView(view, newDataListeners);
        });
      }
      prevDataListenersRef.current = dataListeners;
    };

    const recreateView = () => {
      // console.log('recreating view');
      clearView();
      createView();
    };

    useLayoutEffect(() => {
      // Create the view if it doesn't exist yet
      if (state.viewPromise === undefined) {
        // console.log('Recreating view because promise is undefined');
        recreateView();
      }

      // Only create a new view if necessary
      const specChanges =
        state.spec === null ? false : computeSpecChanges(spec, state.spec);
      // console.log('spec changes', specChanges);

      // Wrap this in modifyView since we only want to recreate the view
      // if the view already exists
      modifyView(() => {
        if (specChanges && specChanges.isExpensive) {
          recreateView();
        }
      });

      // Specify how to clean up after this effect:
      // This function is called when the component unmounts
      return function cleanup() {
        // console.log('CLEANUP: from layout', state);
        clearView();
      };
    }, [spec]);

    useEffect(() => {
      const prevOptions = prevOptionsRef.current;
      const fieldSet = getUniqueFieldNames([options, prevOptions]);

      if (Array.from(fieldSet).some((f) => options[f] !== prevOptions[f])) {
        recreateView();
      }
      prevOptionsRef.current = options;
    }, [options]);

    // Listen to dataset changes and update data
    useEffect(() => {
      modifyView((view) => {
        console.debug('VEGA NEW DATA');

        let prevData = prevDataRef.current;
        let recreate = false;
        if (view && data) {
          Object.keys(data).forEach((name) => {
            if (
              (!Object.prototype.hasOwnProperty.call(prevData, name) ||
                data[name] !== prevData[name]) &&
              recreateOnDatasets.includes(name)
            ) {
              recreate = true;
            }
          });
        }

        if (recreate) {
          recreateView();
        } else {
          updateData(view);
        }
      });
    }, [data]);

    //Update Vega dimensions without redrawing the whole thing
    useEffect(() => {
      modifyView((view) => {
        // console.log('change width/height', width, height);
        if (width !== undefined) {
          view.width(width);
        }
        if (height !== undefined) {
          view.height(height);
        }
        view.run();
      });
    }, [width, height]);

    // Listen to changes in signals passed via. props
    useEffect(() => {
      // console.log('Passed signals:', signals);
      modifyView((view) => {
        let changed = false;
        let doUpdateData = false;
        Object.keys(signals).forEach((signalName) => {
          let currentSignalVal = view.signal(signalName);
          // console.log('current', currentSignalVal, 'new', signals[signalName]);
          // Only update if the signal is different
          if (!equal(currentSignalVal, signals[signalName])) {
            view.signal(signalName, signals[signalName]);
            changed = !cheapSignals.includes(signalName);
            doUpdateData = updateDataSignals.includes(signalName);
          }
        });
        if (changed) {
          // console.log('Signals changed, re-running view...');
          if (doUpdateData) {
            // console.log('Updating data as well');
            prevDataRef.current = {};
            updateData(view);
          } else {
            view.resize().run();
          }
        }
      });
    }, [signals]);

    useEffect(() => {
      updateSignalListeners();
    }, [signalListeners]);

    useEffect(() => {
      updateDataListeners();
    }, [dataListeners]);

    const downloadBlobURL = (url, filename) => {
      var link = document.createElement('a');
      link.setAttribute('href', url);
      link.setAttribute('target', '_blank');
      link.setAttribute('download', filename);
      link.dispatchEvent(new MouseEvent('click'));
    };

    // Download handlers via. refs
    useImperativeHandle(ref, () => ({
      runWhenComplete: (callback) => {
        modifyView((view) => {
          view.runAfter(callback);
        });
      },
      getData: (name, callback) => {
        modifyView((view) => {
          callback(view.data(name));
        });
      },
      downloadImage: (type, filename, scaleFactor = 1) => {
        // type = 'png' or 'svg'
        // scaleFactor = 1 (default)

        modifyView((view) => {
          view
            .toImageURL(type, scaleFactor)
            .then((url) => {
              downloadBlobURL(url, filename);
            })
            .catch((error) => console.error(error));
        });
      },
    }));

    return <div ref={containerRef} className={className} style={style} />;
  }
);

VegaEmbed.propTypes = {
  mode: PropTypes.oneOf(['vega', 'vega-lite']),
  height: PropTypes.number,
  width: PropTypes.number,

  className: PropTypes.string,
  spec: PropTypes.object.isRequired,
  data: PropTypes.object,
  recreateOnDatasets: PropTypes.arrayOf(PropTypes.string),
  signals: PropTypes.object,
  cheapSignals: PropTypes.arrayOf(PropTypes.string),
  updateDataSignals: PropTypes.arrayOf(PropTypes.string),
  signalListeners: PropTypes.object, // key -> value (function)
  dataListeners: PropTypes.object, // key -> value (function)
  style: PropTypes.object,
  onNewView: PropTypes.oneOfType([PropTypes.func, PropTypes.oneOf([null])]),
  onComplete: PropTypes.oneOfType([PropTypes.func, PropTypes.oneOf([null])]),
  onError: PropTypes.func,
};
VegaEmbed.defaultProps = {
  mode: 'vega',
  className: 'vega-embed',
  data: {},
  recreateOnDatasets: [],
  signals: {},
  cheapSignals: [],
  updateDataSignals: [],
  signalListeners: {},
  dataListeners: {},
  style: {},
  onNewView: null,
  onComplete: null,
  onError: (error) => {
    console.error(error);
    return null;
  },
};

export default VegaEmbed;
