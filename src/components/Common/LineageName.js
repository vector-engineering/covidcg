import React from 'react';
import PropTypes from 'prop-types';
import { Name } from './LineageName.styles';

const determineWHOLabel = (name) => {
  const regExpressions = {
    '^B.1.1.7.?': 'Alpha',
    '^Q.?': 'Alpha',
    '^B.1.351.?': 'Beta',
    '^B.1.1.281.?': 'Gamma',
    '^P.1.?': 'Gamma',
    '^B.1.617.2.?': 'Delta',
    '^AY.?': 'Delta',
    '^B.1.1.529.?': 'Omicron',
    '^BA.1.?': 'Omicron',
    '^C.37.?': 'Lambda',
    '^B.1.621.?': 'Mu',
    '^BB.2.?': 'Mu',
  };

  for (let key of Object.keys(regExpressions)) {
    const regEx = new RegExp(key, 'i');
    if (name.match(regEx)) {
      return regExpressions[key];
    }
  }

  return null;
};

export const LineageName = ({ name, selected }) => {
  const label = determineWHOLabel(name);
  return (
    <Name selected={selected} whoLabel={label}>
      {name}
    </Name>
  );
};
LineageName.propTypes = {
  name: PropTypes.string.isRequired,
  selected: PropTypes.bool,
};
LineageName.defaultProps = {
  selected: false,
};
