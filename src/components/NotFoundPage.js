import React from 'react';
import PropTypes from 'prop-types';
import { Link } from 'mobx-router';
import routes from '../routes';
import { connect } from '../stores/connect';

const NotFoundPage = (props) => {
  return (
    <div>
      <h4>404 Page Not Found</h4>
      <Link router={props.router} route={routes.home}>
        Go back to homepage
      </Link>
    </div>
  );
};

NotFoundPage.propTypes = {
  router: PropTypes.object.isRequired,
};

export default connect(NotFoundPage);
