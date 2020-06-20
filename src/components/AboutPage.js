import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { Link } from 'mobx-router';
import routes from '../routes';
import { connect } from '../stores/connect';

const AboutPageContainer = styled.div`
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  justify-content: flex-start;
`;
const Header = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  margin-bottom: 10px;
  padding: 20px;
  border-bottom: 1px solid #aaa;
  width: 100%;
`;
const Title = styled.h1`
  margin: 0px;
  margin-right: 20px;
  font-weight: 700;
  font-size: 2em;
  border-right: 1px solid #aaa;
  padding-right: 20px;
`;
const Content = styled.div`
  max-width: 1000px;
  padding: 20px;
`;
const ImageRow = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-start;
  a {
    margin-right: 30px;
  }
`;

// Since this component is simple and static, there's no parent container for it.
const AboutPage = (props) => {
  return (
    <AboutPageContainer>
      <Header>
        <Title>About</Title>
        <Link route={routes.home} router={props.router}>
          Back to Home
        </Link>
      </Header>

      <Content>
        <p>
          This app was made by....
          <br />
          <ImageRow>
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
          </ImageRow>
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
        <ImageRow>
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
        </ImageRow>
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
        https://github.com/matiassingers/emoji-flags
      </Content>
    </AboutPageContainer>
  );
};

AboutPage.propTypes = {
  router: PropTypes.object.isRequired,
};

export default connect(AboutPage);
