import React from 'react';
import PropTypes from 'prop-types';

const ExternalLink = ({ children, href }) => {
  return (
    <a href={href} target="_blank" rel="noopener noreferrer">
      {children}
    </a>
  );
};
ExternalLink.propTypes = {
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node,
  ]),
  href: PropTypes.string.isRequired,
};
ExternalLink.defaultProps = {
  children: '',
};

export default ExternalLink;
