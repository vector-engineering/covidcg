import React from 'react';
import { Link } from 'mobx-router';
import routes from '../routes';
import { connect } from '../stores/connect';

const NotFoundPage = (props) => {
  return (
    <div>
      <h4>404 Page Not Found</h4>
      <Link router={props.router} view={routes.home}>
        Go back to homepage
      </Link>
    </div>
  );
};

export default connect(NotFoundPage);
