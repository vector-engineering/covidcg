import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import Modal from 'react-modal';

import CGLogo from '../../assets/images/cg_logo_v13.png';

import ExternalLink from '../Common/ExternalLink';

const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;
  justify-content: flex-start;

  padding: 0px 10px;
`;
const Header = styled.div`
  display: flex;
  flex-direction: row;
  align-items: flex-start;

  margin-bottom: 10px;
  padding-top: 10px;

  .title {
    h2 {
      margin-bottom: 0px;
      margin-top: 0px;
      margin-left: 5px;
    }
  }
  .spacer {
    flex-grow: 1;
  }
  .close-button {
  }
`;
const Content = styled.div`
  font-size: 1em;
  font-weight: normal;
`;

Modal.setAppElement('#app');
const NOOP = () => {};

const SplashScreenContent = ({ onRequestClose }) => {
  return (
    <Wrapper>
      <Header>
        <div className="title">
          <img height={75} src={CGLogo}></img>
          <h2>COVID-19 CoV Genetics</h2>
        </div>
        <div className="spacer"></div>
        <button className="close-button" onClick={onRequestClose}>
          Dismiss
        </button>
      </Header>
      <Content>
        <p>
          <b>
            The COVID-19 CoV Genetics browser was designed to empower diverse
            projects on SARS-CoV-2 transmission, evolution, emergence, immune
            interactions, diagnostics, therapeutics, vaccines, and tracking of
            interventions.
          </b>
        </p>
        <p>
          Being able to track SARS-CoV-2 evolution - updated daily using{' '}
          <ExternalLink href="https://www.gisaid.org/">GISAID</ExternalLink>{' '}
          virus data! - according to user-defined locations, dates, and
          mutations is very important for scientists and public health
          professionals, as well as developers of vaccines, diagnostics, and
          therapeutics.
        </p>
        <p>
          COVID-19 CG helps users to quickly find answers to questions including
          but not limited to:
        </p>
        <ol>
          <li>
            Which SARS-CoV-2 lineages are present in a given city or region
            within a user-specified period of time?
          </li>
          <li>
            Which SARS-CoV-2 variants should I test my therapeutic, antibody, or
            diagnostic on before implementation in a specific region?
          </li>
          <li>
            What are the community outcomes after particular policies, vaccines,
            or therapeutics are applied in that population?
          </li>
          <li>
            Are there data from transient mutations that can elucidate common
            mechanisms of resistance to acquired immunity? Can this be leveraged
            for vaccine, antibody, or small molecule drug design?
          </li>
        </ol>
        <p>
          Users can view the comprehensive nucleotide and amino acid residue
          variation in their selection to inform their research hypothesis
          generation or anti-COVID-19 product development. For example, COVID-19
          CG enables users to evaluate commonly used or custom primers/probes or
          targets/epitopes based on their location and dates of interest. As
          sequencing efforts start to include details about SARS-CoV-2 isolates,
          users can also sort virus data according to patient characteristics
          such as age, gender, clinical status, isolate type, as well as
          passaging, sequencing, and assembly method.
        </p>
        <p>
          To accelerate COVID-19 research and public health efforts, COVID-19 CG
          will be continually upgraded with new and improved features so that
          users can quickly and reliably pinpoint critical mutations as the
          virus evolves throughout the pandemic.
        </p>
        <p>
          Towards this goal, we strongly advocate that countries continue to
          sequence SARS-CoV-2 isolates from their patients (and infected
          animals) and share this data in a timely manner via GISAID so that
          scientists worldwide are maximally informed about developments in the
          spread of SARS-CoV-2.
        </p>
        <p>
          Reach out to us{' '}
          <ExternalLink href="https://twitter.com/covidcg">
            @covidcg
          </ExternalLink>{' '}
          on twitter.
        </p>
        <p>
          Find out more about how you can use COVID-19 CG in our preprint: URL
          and maybe youtube video
        </p>
      </Content>
    </Wrapper>
  );
};
SplashScreenContent.propTypes = {
  onRequestClose: PropTypes.func,
};
SplashScreenContent.defaultProps = {
  onRequestClose: NOOP,
};

const SplashScreenModal = ({ isOpen, onAfterOpen, onRequestClose }) => {
  return (
    <Modal
      isOpen={isOpen}
      onAfterOpen={onAfterOpen}
      onRequestClose={onRequestClose}
      style={{
        overlay: {
          zIndex: 2,
        },
        content: {
          top: '50%',
          left: '50%',
          right: 'auto',
          bottom: 'auto',
          marginRight: '-50%',
          transform: 'translate(-50%, -50%)',
          maxWidth: '800px',
          zIndex: 3,
        },
      }}
      contentLabel="Splash Screen Modal"
    >
      <SplashScreenContent onRequestClose={onRequestClose} />
    </Modal>
  );
};

SplashScreenModal.propTypes = {
  isOpen: PropTypes.bool.isRequired,
  onAfterOpen: PropTypes.func,
  onRequestClose: PropTypes.func.isRequired,
};

SplashScreenModal.defaultProps = {
  onAfterOpen: NOOP,
};

export default SplashScreenModal;
