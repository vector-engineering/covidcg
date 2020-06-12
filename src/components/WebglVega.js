import React, { useLayoutEffect, useRef } from 'react';
import PropTypes from 'prop-types';

import { toJS } from 'mobx';

const WebglVega = ({ data, spec }) => {
  const target = useRef();

  spec.data.values = toJS(data).slice(0, 10);
  spec.$schema = 'https://vega.github.io/schema/vega-lite/v3.json';

  console.log('ourspec: ', JSON.stringify(spec));

  useLayoutEffect(() => {
    //console.log('layout effect');
    window.d3.json(
      'https://vega.github.io/schema/vega-lite/v3.json',
      (schema) => {
        console.log('hello: ', spec, schema);
        new window.vega.View(window.vega.parse(spec))
          .initialize(target.current)
          .renderer('webgl')
          .run();
      }
    );

    console.log(target);
  }, []);
  return <div ref={target} style={{ width: '100%', height: '800px' }}></div>;
};

WebglVega.propTypes = {
  data: PropTypes.arrayOf(PropTypes.object).isRequired,
  spec: PropTypes.object.isRequired,
};

export default WebglVega;
