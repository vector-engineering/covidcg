import React from 'react';
import styled from 'styled-components';

import ExternalLink from '../Common/ExternalLink';

// import ReactSlingshotImage from '../../assets/images/react_slingshot.png';
// import ReactLogo from '../../assets/images/React-icon.svg';
// import MobXLogo from '../../assets/images/mobx.png';
// import IDLLogo from '../../assets/images/idl-logo.png';
// import UBCLogo from '../../assets/images/ubc_logo.png';
import VELogo from '../../assets/images/VE_logo_new.png';
import BroadLogo from '../../assets/images/BroadLogo_RGB_forDigital.png';
import GISAIDLogo from '../../assets/images/gisaid_logo.png';

import {
  TabContainer,
  Content,
  ContentSection,
  ImageRow,
} from './TextTab.styles';

const VELogoText = styled.span`
  font-size: 1.5em;
  font-weight: bold;
`;

// Since this component is simple and static, there's no parent container for it.
const AboutTab = () => {
  return (
    <TabContainer>
      <Content>
        {/*
        <TOC>
          <span className="toc-title">Table of Contents</span>
          <ul className="toc-list">
            <li className="toc-item">
              <a title="contributors" href="#contributors">
                Contributors
              </a>
            </li>
          </ul>
        </TOC>*/}
        <ContentSection>
          <a id="contributors"></a>
          <span className="section-title">COVID-CG is developed by</span>

          <div className="content-text">
            <p>
              COVID CoV Genomics (CG) was developed in the{' '}
              <ExternalLink href="https://vector-engineering.github.io/">
                Vector Engineering Lab
              </ExternalLink>{' '}
              (
              <ExternalLink href="https://www.broadinstitute.org/stanley">
                Stanley Center for Psychiatric Research
              </ExternalLink>
              ,{' '}
              <ExternalLink href="https://www.broadinstitute.org/">
                Broad Institute of MIT and Harvard
              </ExternalLink>
              ), by:
            </p>

            <p>
              <b>Albert Tian Chen</b> is an associate computational biologist in
              the Deverman group who designed, built, and will continue to
              expand the COVID-19 CoV Genetics browser. Albert has previous
              experience building websites for startups and interactive data
              visualization tools for research groups.
            </p>
            <p>
              <b>Kevin Altschuler</b> is an experienced web developer with
              previous positions at Google and Uber Advanced Technologies Group
              (ATG). Kevin has architected the application to streamline
              development and improve performance. (
              <ExternalLink href="https://www.linkedin.com/in/kevinaltschuler/">
                LinkedIn Profile
              </ExternalLink>
              )
            </p>
            <p>
              <b>Dr. Alina Chan</b> is a postdoc in the Deverman group who asked
              too many questions about SARS-CoV-2 and realized that we needed a
              browser like COVID-19 CG.
            </p>
            <p>
              <b>Dr. Shing Hei Zhan</b> is a recent graduate from the University
              of British Columbia’s Department of Zoology &amp; the Biodiversity
              Research Centre. Shing is also co-founder and lead bioinformatics
              scientist at Fusion Genomics Corporation, which develops molecular
              diagnostic assays for infectious diseases, including COVID-19.
            </p>
            <p>
              <b>Dr. Benjamin Deverman</b> is the director of the vector
              engineering research group at the Stanley Center for Psychiatric
              Research at the Broad Institute of MIT and Harvard. The Deverman
              Lab develops innovative gene delivery solutions to expand the
              impact of gene therapy and for studying the central nervous
              system, with the aim of uncovering new avenues for treating
              psychiatric disorders.
            </p>

            <p>
              Contact the authors by email:{' '}
              <a href="mailto:covidcg@broadinstitute.org">
                covidcg@broadinstitute.org
              </a>
            </p>

            <ImageRow>
              <ExternalLink href="https://vector.engineering" showIcon={false}>
                <img src={VELogo} height="40"></img>
              </ExternalLink>
              <VELogoText>Vector Engineering Lab</VELogoText>
              {/* <ExternalLink href="https://www.broadinstitute.org/">
                <img
                  src={BroadLogo}
                  height="40"
                  style={{ marginLeft: '20px' }}
                ></img>
              </ExternalLink> */}
            </ImageRow>
            <ImageRow>
              {/* <ExternalLink href="https://www.broadinstitute.org/stanley">
                  <img
                    src="https://storage.googleapis.com/ve-public/covid_ui/assets/img/StanleyCenterLogo_RGB_forDigital.png"
                    height="50"
                  ></img>
                </ExternalLink> */}
            </ImageRow>
            {/* <ImageRow>
                <ExternalLink href="https://www.ubc.ca/">
                  <img src={UBCLogo} height="60"></img>
                </ExternalLink>
              </ImageRow> */}
          </div>
        </ContentSection>

        <ContentSection>
          <a id="sequence-data"></a>
          <span className="section-title">Data enabling COVID-CG</span>

          <div className="content-text">
            <p>
              We are extremely grateful to the{' '}
              <ExternalLink href="https://www.gisaid.org/">
                GISAID Initiative
              </ExternalLink>{' '}
              and all its data contributors, i.e. the Authors from the
              Originating laboratories responsible for obtaining the specimens
              and the Submitting laboratories where genetic sequence data were
              generated and shared via the GISAID Initiative, on which this
              research is based.
            </p>
            {/* <ImageRow>
              <ExternalLink href="https://gisaid.org">
                <img src={GISAIDLogo} height="60" />
              </ExternalLink>
            </ImageRow> */}
            {/* <iframe
              id="gisaid-pub-pdf"
              style={{ height: '30em' }}
              width="100%"
              src={
                'https://ve-public.storage.googleapis.com/eurosurv-22-30494-1.pdf#page=1&zoom=100'
              }
            /> */}
            <p>
              Elbe, S., and Buckland-Merrett, G. (2017) Data, disease and
              diplomacy: GISAID’s innovative contribution to global health.{' '}
              <i>Global Challenges</i>, 1:33-46. DOI:
              <ExternalLink
                href="https://doi.org/10.1002/gch2.1018"
                showIcon={false}
              >
                10.1002/gch2.1018
              </ExternalLink>{' '}
              PMCID:{' '}
              <ExternalLink
                href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6607375/"
                showIcon={false}
              >
                31565258
              </ExternalLink>
            </p>
          </div>
        </ContentSection>

        <ContentSection>
          <a id="citing-covid-cg"></a>
          <span className="section-title">Citing COVID-19 CG:</span>
          <div className="content-text">
            {/* <p>
                Chen AT, Altschuler K, Chan AY, Zhan SH, Deverman BE (2020).
                COVID-19 CG: Tracking SARS-CoV-2 by mutation, location, and date
                of interest. <i>bioRxiv</i>. DOI: ...
              </p> */}
            <p>
              Users are encouraged to share, download, and further analyze data
              from this site. Plots can be downloaded as PNG or SVG files, and
              the data powering the plots and tables can be downloaded as well.
              Please attribute any data/images to{' '}
              <a href="https://covidcg.org">covidcg.org</a>.
            </p>
            <p>
              Note: When using results from these analyses in your manuscript,
              ensure that you acknowledge the contributors of data, i.e.{' '}
              <i>
                We gratefully acknowledge all the Authors from the Originating
                laboratories responsible for obtaining the speciments and the
                Submitting laboratories where genetic sequence data were
                generated and shared via the GISAID Initiative, on which this
                research is based.
              </i>
            </p>
            <p>and cite the following reference(s):</p>
            <p>
              Shu, Y., McCauley, J. (2017) GISAID: Global initiative on sharing
              all influenza data – from vision to reality.{' '}
              <i>EuroSurveillance</i>, 22(13) DOI:
              <ExternalLink
                href="https://doi.org/10.2807/1560-7917.ES.2017.22.13.30494"
                showIcon={false}
              >
                10.2807/1560-7917.ES.2017.22.13.30494
              </ExternalLink>{' '}
              PMCID:{' '}
              <ExternalLink
                href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5388101/"
                showIcon={false}
              >
                PMC5388101
              </ExternalLink>
            </p>
          </div>
        </ContentSection>

        {/* <ContentSection>
          <a id="reference-data"></a>
          <span className="section-title">Reference Sequence</span>

          <div className="content-block">
            <div className="content-text">
              <p>
                All hCoV-19 sequences are aligned to the WIV04 reference
                sequence (GISAID: EPI_ISL_402124). Gene and protein ORFs were taken from the GenBank page.
              </p>
            </div>
            <div className="content-images"></div>
          </div>
        </ContentSection> */}
        {/*
        <ContentSection>
          <a id="primer-data"></a>
          <span className="section-title">Primer / Probe Data</span>

          <div className="content-block">
            <div className="content-text"></div>
            <div className="content-images"></div>
          </div>
        </ContentSection>
        */}

        {/* <ContentSection>
          <a id="code"></a>
          <span className="section-title">Open-source Code</span>

          <p>
            This project is built off of the following open-source projects and
            its contributors
          </p>

          <div className="content-block">
            <div className="content-text">
              <p>
                This app is built on the{' '}
                <ExternalLink href="https://reactjs.org/">
                  React.js
                </ExternalLink>{' '}
                framework, and was initially made from the{' '}
                <ExternalLink href="https://github.com/coryhouse/react-slingshot">
                  React-Slingshot starter kit
                </ExternalLink>
                .{' '}
                <ExternalLink href="https://mobx.js.org/README.html">
                  MobX
                </ExternalLink>{' '}
                is used for internal state management.
              </p>
            </div>
            <div className="content-images">
              <ImageRow>
                <ExternalLink href="https://reactjs.org/">
                  <img src={ReactLogo} height="50"></img>
                </ExternalLink>
                <ExternalLink href="https://github.com/coryhouse/react-slingshot">
                  <img src={ReactSlingshotImage} height="40"></img>
                </ExternalLink>
                <ExternalLink href="https://mobx.js.org/README.html">
                  <img src={MobXLogo} height="50"></img>
                </ExternalLink>
              </ImageRow>
            </div>
          </div>

          <div className="content-block">
            <div className="content-text">
              <p>
                This project uses{' '}
                <ExternalLink href="https://vega.github.io/">Vega</ExternalLink>{' '}
                from the{' '}
                <ExternalLink href="https://idl.cs.washington.edu/">
                  University of Washington Interactive Data Lab
                </ExternalLink>{' '}
                to generate deep and interactive data visualizations.{' '}
                <ExternalLink href="https://github.com/vega/react-vega">
                  react-vega
                </ExternalLink>{' '}
                is used to interface Vega with our UI.
              </p>
            </div>
            <div className="content-images">
              <ExternalLink href="https://idl.cs.washington.edu/">
                <img src={IDLLogo} height="50"></img>
              </ExternalLink>
            </div>
          </div>

          <div className="content-block">
            <div className="content-text">
              <p>
                Tree-selectors are built from{' '}
                <ExternalLink href="https://github.com/dowjones/react-dropdown-tree-select">
                  dowjones/react-dropdown-tree-select
                </ExternalLink>
                . Flags in the location selector are taken from{' '}
                <ExternalLink href="https://github.com/matiassingers/emoji-flags">
                  matiassingers/emoji-flags
                </ExternalLink>
              </p>
            </div>
            <div className="content-images"></div>
          </div>

          <div className="content-block">
            <div className="content-text">
              <p>
                Tables are rendered with{' '}
                <ExternalLink href="https://adazzle.github.io/react-data-grid/">
                  adazzle/react-data-grid
                </ExternalLink>
              </p>
            </div>
            <div className="content-images"></div>
          </div>

          <div className="content-block">
            <div className="content-text"></div>
            <div className="content-images"></div>
          </div>
        </ContentSection> */}
      </Content>
    </TabContainer>
  );
};

export default AboutTab;
