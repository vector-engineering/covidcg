import React from 'react';
import PropTypes from 'prop-types';
import '../styles/about-page.scss';
import { Link } from 'mobx-router';
import routes from '../routes';
import { connect } from '../stores/connect';

// Since this component is simple and static, there's no parent container for it.
const AboutPage = (props) => {
  return (
    <div className="about-page">
      <div className="header">
        <h1>About</h1>
        <Link route={routes.home} router={props.router}>
          Back to Home
        </Link>
      </div>

      <div className="content">
        <p>
          This app was made by....
          <br />
          <div className="image-row">
            <a
              href="https://www.broadinstitute.org/stanley-center-psychiatric-research/vector-engineering"
              target="_blank"
              rel="noopener noreferrer"
            >
              <img
                src="https://storage.googleapis.com/ve-public/covid_ui/assets/img/ve_logo.png"
                height="60"
              ></img>
            </a>
            <a
              href="https://www.broadinstitute.org/stanley"
              target="_blank"
              rel="noopener noreferrer"
            >
              <img
                src="https://storage.googleapis.com/ve-public/covid_ui/assets/img/StanleyCenterLogo_RGB_forDigital.png"
                height="60"
              ></img>
            </a>
            <a
              href="https://www.broadinstitute.org/"
              target="_blank"
              rel="noopener noreferrer"
            >
              <img
                src="https://storage.googleapis.com/ve-public/covid_ui/assets/img/BroadLogo_RGB_forDigital.png"
                height="60"
              ></img>
            </a>
          </div>
          <a
            href="https://www.broadinstitute.org/stanley-center-psychiatric-research/vector-engineering"
            target="_blank"
            rel="noopener noreferrer"
          >
            Vector Engineering Group
          </a>
          ,{' '}
          <a
            href="https://www.broadinstitute.org/stanley"
            target="_blank"
            rel="noopener noreferrer"
          >
            Stanley Center for Psychiatric Research
          </a>
          ,{' '}
          <a
            href="https://www.broadinstitute.org/"
            target="_blank"
            rel="noopener noreferrer"
          >
            Broad Institute of MIT and Harvard
          </a>
        </p>

        <h2>Attributions</h2>

        {/* https://github.com/dowjones/react-dropdown-tree-select
        https://github.com/vega/react-vega
        underscore.js
        VEGA
        react-data-table-component https://github.com/jbetancur/react-data-table-component#readme */}

        <h3>Data</h3>

        <div className="image-row">
          <a
            href="https://www.gisaid.org/"
            target="_blank"
            rel="noopener noreferrer"
          >
            <img
              src="https://storage.googleapis.com/ve-public/covid_ui/assets/img/gisaid.png"
              height="60"
            ></img>
          </a>
        </div>
        <p>
          <a
            href="https://www.gisaid.org/"
            target="_blank"
            rel="noopener noreferrer"
          >
            GISAID
          </a>{' '}
          for all the sequences
        </p>

        <h3>Code</h3>

        <p>
          This app was made from the{' '}
          <a
            href="https://github.com/coryhouse/react-slingshot"
            target="_blank"
            rel="noopener noreferrer"
          >
            React-Slingshot starter kit
          </a>
          .
        </p>
      </div>
    </div>
  );
};

AboutPage.propTypes = {
  router: PropTypes.object.isRequired,
};

export default connect(AboutPage);
