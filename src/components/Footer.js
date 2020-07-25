import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { NOOP } from '../constants/functions';

const FooterContainer = styled.div`
  margin-top: auto;
  display: flex;
  background-color: #f8f8f8;

  margin-left: -10px;
  padding: 5px;
  padding-left: 20px;
  border-top: 1px solid #ccc;

  font-size: 0.85rem;

  .gisaid-daa {
    margin-right: 10px;
    padding-right: 10px;
    border-right: 1px solid #aaa;
  }
`;

const Footer = ({ openModal }) => {
  return (
    <FooterContainer>
      <div className="gisaid-daa">
        Data use subject to the{' '}
        <a
          href="https://www.gisaid.org/"
          target="_blank"
          rel="noopener noreferrer"
        >
          GISAID
        </a>{' '}
        EpiCov™{' '}
        <a
          href="https://www.gisaid.org/registration/terms-of-use/"
          target="_blank"
          rel="noopener noreferrer"
        >
          Database Access Agreement
        </a>
      </div>
      <a href="#" onClick={openModal}>
        Show Splash Screen
      </a>
    </FooterContainer>
  );
};
Footer.propTypes = {
  openModal: PropTypes.func,
};
Footer.defaultProps = {
  openModal: NOOP,
};

export default Footer;
