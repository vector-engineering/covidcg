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

const VegaEmbed = ({
  className,
  spec,
  data,
  signalListeners,
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

  const prevOptionsRef = useRef(null);
  const prevSpecRef = useRef(null);
  const prevWidthRef = useRef(null);
  const prevHeightRef = useRef(null);
  const prevSignalListenersRef = useRef(null);

  const handleError = (error) => {
    onError(error);
    // eslint-disable-next-line no-console
    console.warn(error);
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

  const updateSingleDatasetInView = (view, name, value) => {
    if (value) {
      if (isFunction(value)) {
        value(view.data(name));
      } else {
        view.change(
          name,
          vega
            .changeset()
            .remove(() => true)
            .insert(value)
        );
      }
    }
  };

  const updateMultipleDatasetsInView = (view, data) => {
    Object.keys(data).forEach((name) => {
      updateSingleDatasetInView(view, name, data[name]);
    });
  };

  const createView = () => {
    // console.log('CREATE VIEW');
    if (containerRef.current) {
      const finalSpec = combineSpecWithDimension({ spec, width, height });
      viewPromise = vegaEmbed(containerRef.current, finalSpec, options)
        .then(({ view }) => {
          if (addSignalListenersToView(view, signalListeners)) {
            view.run();
          }
          return view;
        })
        .catch(handleError);

      if (onNewView) {
        modifyView(onNewView);
      }
    }
  };

  const modifyView = (action) => {
    if (viewPromise) {
      viewPromise
        .then((view) => {
          if (view) {
            action(view);
          }

          return true;
        })
        .catch(handleError);
    }
  };

  const clearView = () => {
    modifyView((view) => {
      view.finalize();
    });
    viewPromise = undefined;

    return this;
  };

  useEffect(() => {
    prevSpecRef.current = spec;
    prevOptionsRef.current = options;
    prevWidthRef.current = width;
    prevHeightRef.current = height;
    prevSignalListenersRef.current = signalListeners;
  });

  useEffect(() => {
    // console.log('new data');
    // console.log(data);

    if (data && Object.keys(data).length > 0) {
      modifyView((view) => {
        updateMultipleDatasetsInView(view, data);
        view.resize().run();
      });
    }
  }, [data]);

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

    fieldSet.delete('className');
    fieldSet.delete('signalListeners');
    fieldSet.delete('spec');
    fieldSet.delete('style');
    fieldSet.delete('width');
    fieldSet.delete('height');

    const prevSpec = prevSpecRef.current;
    const prevWidth = prevWidthRef.current;
    const prevHeight = prevHeightRef.current;

    // Only create a new view if necessary
    if (Array.from(fieldSet).some((f) => options[f] !== prevOptions[f])) {
      clearView();
      createView();
    } else {
      // BUG: width/height changes only are being marked as expensive, for some reason
      const specChanges = computeSpecChanges(
        combineSpecWithDimension({ spec, width, height }),
        combineSpecWithDimension({
          spec: prevSpec,
          width: prevWidth,
          height: prevHeight,
        })
      );

      const newSignalListeners = signalListeners;
      const oldSignalListeners = prevSignalListenersRef.current;

      if (specChanges) {
        if (specChanges.isExpensive) {
          clearView();
          createView();
        } else {
          const areSignalListenersChanged = !shallowEqual(
            newSignalListeners,
            oldSignalListeners
          );
          modifyView((view) => {
            if (specChanges.width !== false) {
              view.width(specChanges.width);
            }
            if (specChanges.height !== false) {
              view.height(specChanges.height);
            }
            if (areSignalListenersChanged) {
              removeSignalListenersFromView(view, oldSignalListeners);
              addSignalListenersToView(view, newSignalListeners);
            }

            view.run();
          });
        }
      } else if (!shallowEqual(newSignalListeners, oldSignalListeners)) {
        modifyView((view) => {
          removeSignalListenersFromView(view, oldSignalListeners);
          addSignalListenersToView(view, newSignalListeners);
          view.run();
        });
      }
    }

    // Specify how to clean up after this effect:
    // This function is called when the component unmounts
    return function cleanup() {
      clearView();
    };
  });

  return <div ref={containerRef} className={className} style={style} />;
};

VegaEmbed.propTypes = {
  mode: PropTypes.oneOf(['vega', 'vega-lite']),
  height: PropTypes.number,
  width: PropTypes.number,

  className: PropTypes.string,
  spec: PropTypes.object.isRequired,
  data: PropTypes.object.isRequired,
  signalListeners: PropTypes.object, // key -> value (function)
  style: PropTypes.object,
  onNewView: PropTypes.func,
  onError: PropTypes.func,
};
VegaEmbed.defaultProps = {
  mode: 'vega',
  className: 'vega-embed',
  signalListeners: {},
  style: {},
  onNewView: NOOP,
  onError: (error) => {
    console.error(error);
    return null;
  },
};

export default VegaEmbed;
