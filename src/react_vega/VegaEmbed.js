import React, {
  createRef,
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

const VegaEmbed = forwardRef(
  (
    {
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
    },
    ref
  ) => {
    const containerRef = createRef();

    const [state, setState] = useState({
      initialized: false,
      viewPromise: undefined,
    });

    const prevOptionsRef = useRef({});
    const prevSpecRef = useRef({});
    const prevDataRef = useRef({});
    const prevSignalListenersRef = useRef({});
    const prevDataListenersRef = useRef({});

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
        // vega-embed options: https://github.com/vega/vega-embed#options
        let viewPromise = vegaEmbed(containerRef.current, finalSpec, options)
          .then(({ view }) => {
            // Unset previous data to force push the data into the new view
            prevDataRef.current = {};
            updateData(view);

            return { view };
          })
          .then(({ view }) => {
            addSignalListenersToView(view, signalListeners);
            addDataListenersToView(view, dataListeners);
            view.height(3); // TODO: I have no idea why this works or what it does
            view.run();
            return { view };
          })
          .catch(handleError);

        setState({ ...state, viewPromise });

        if (onNewView !== null) {
          modifyView(onNewView);
        }
      }
    };

    const modifyView = (action) => {
      // console.log('MODIFY VIEW');
      if (state.viewPromise !== undefined) {
        state.viewPromise
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
      setState({ ...state, viewPromise: undefined });

      prevOptionsRef.current = {};
      prevSpecRef.current = {};
      prevDataRef.current = {};
      prevSignalListenersRef.current = {};
      prevDataListenersRef.current = {};

      return this;
    };

    const updateData = (view) => {
      if (view && data) {
        let prevData = prevDataRef.current;

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
            // console.log('Changing datasets', name);
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

        // let promise = view.runAsync();
        // promise.then((view) => {
        //   console.log('finished');
        //   console.log(view.getState());
        // });
        view.runAsync();
      }
      prevDataRef.current = data;
    };

    // Listen to changes in signals passed via. props
    useEffect(() => {
      // console.log('Passed signals:', signals);
      modifyView((view) => {
        Object.keys(signals).forEach((signalName) => {
          let currentSignalVal = view.signal(signalName);
          // console.log('current', currentSignalVal, 'new', signals[signalName]);
          // Only update if the signal is different
          if (currentSignalVal != signals[signalName]) {
            view.signal(signalName, { group: signals[signalName] });
          }
        });
      });
    }, [signals]);

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
    };

    const recreateView = () => {
      // If the view is currently being re-rendered, then skip this...
      if (state.viewPromise !== undefined) {
        // console.log('in the middle of a re-render');
        // return;
      }

      // console.log('recreating view');
      clearView();
      createView();
      // updateSignalListeners();
      // updateDataListeners();
      // updateData();
    };

    useLayoutEffect(() => {
      if (state.viewPromise === undefined) {
        recreateView();
        // updateData(state.view);
      }

      const prevSpec = prevSpecRef.current;

      // Only create a new view if necessary
      const specChanges = computeSpecChanges(spec, prevSpec);
      // console.log(spec, prevSpec);
      // console.log('spec changes', specChanges);

      if (specChanges && specChanges.isExpensive) {
        // console.log('Triggering recreate because of spec changes');
        recreateView();
      }

      // Specify how to clean up after this effect:
      // This function is called when the component unmounts
      return function cleanup() {
        // console.log('CLEANUP: from layout');
        clearView();
      };
    }, [spec]);

    useEffect(() => {
      const prevOptions = prevOptionsRef.current;
      const fieldSet = getUniqueFieldNames([options, prevOptions]);
      // console.log(options, prevOptions, fieldSet);

      if (Array.from(fieldSet).some((f) => options[f] !== prevOptions[f])) {
        // console.log('Trigger recreate because of option changes');
        recreateView();
      }
    }, [options]);

    // Listen to dataset changes and update data
    useEffect(() => {
      //console.log('NEW DATA');
      modifyView((view) => {
        updateData(view);
      });
    }, [data]);

    //Update Vega dimensions without redrawing the whole thing
    useEffect(() => {
      modifyView((view) => {
        console.log('change width/height', width, height);
        if (width !== undefined) {
          view.width(width);
        }
        if (height !== undefined) {
          view.height(height);
        }
        view.run();
      });
    }, [width, height]);

    useEffect(() => {
      updateSignalListeners();
    }, [signalListeners]);

    useEffect(() => {
      updateDataListeners();
    }, [dataListeners]);

    useEffect(() => {
      prevSpecRef.current = spec;
      prevOptionsRef.current = options;
      prevSignalListenersRef.current = signalListeners;
      prevDataListenersRef.current = dataListeners;
    }, [spec, options, signalListeners, dataListeners]);

    const downloadBlobURL = (url, filename) => {
      var link = document.createElement('a');
      link.setAttribute('href', url);
      link.setAttribute('target', '_blank');
      link.setAttribute('download', filename);
      link.dispatchEvent(new MouseEvent('click'));
    };

    // Download handlers via. refs
    useImperativeHandle(ref, () => ({
      downloadImage: (type, filename, scaleFactor = 1) => {
        // console.log('DOWNLOAD PNG');
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
  data: PropTypes.object.isRequired,
  signals: PropTypes.object,
  signalListeners: PropTypes.object, // key -> value (function)
  dataListeners: PropTypes.object, // key -> value (function)
  style: PropTypes.object,
  onNewView: PropTypes.oneOf([PropTypes.func, null]),
  onError: PropTypes.func,
};
VegaEmbed.defaultProps = {
  mode: 'vega',
  className: 'vega-embed',
  signals: {},
  signalListeners: {},
  dataListeners: {},
  style: {},
  onNewView: null,
  onError: (error) => {
    console.error(error);
    return null;
  },
};

export default VegaEmbed;
