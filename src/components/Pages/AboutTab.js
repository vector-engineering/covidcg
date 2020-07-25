import React from 'react';
import styled from 'styled-components';

import ExternalLink from '../ExternalLink';

import ReactSlingshotImage from '../../assets/images/react_slingshot.png';
import ReactLogo from '../../assets/images/React-icon.svg';
import MobXLogo from '../../assets/images/mobx.png';
import IDLLogo from '../../assets/images/idl-logo.png';
import NextstrainLogo from '../../assets/images/nextstrain_logo.png';
import HKUSTLogo from '../../assets/images/HKUST-original_0.svg';
import LANLLogo from '../../assets/images/lanl_logo.svg';
import PangolinLogo from '../../assets/images/pangolin_logo.png';
import UCLLogo from '../../assets/images/university-college-london-ucl-vector-logo.svg';
import COGUKLogo from '../../assets/images/logo-cog-uk.png';
import UBCLogo from '../../assets/images/ubc_logo.png';
import JHULogo from '../../assets/images/jhu_logo.jpg';

const AboutTabContainer = styled.div`
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  justify-content: flex-start;

  background-color: #f8f8f8;
`;

const Content = styled.div`
  flex-grow: 1;
  max-width: 1000px;
  padding: 20px;
  margin: 0 auto;

  background-color: #fff;
`;

// const TOC = styled.div`
//   span.toc-title {
//     font-weight: bold;
//   }
// `;

const ContentSection = styled.div`
  border-top: 1px solid #ccc;
  font-weight: normal;
  font-size: 1em;
  padding-top: 10px;
  padding-bottom: 10px;

  .section-title {
    font-weight: bold;
    font-size: 1.5em;
  }

  .content-block {
    display: flex;
    flex-direction: row;
    align-items: flex-start;

    margin-top: 10px;

    .content-text {
      width: 60%;
      padding-right: 10px;

      .content-subtitle {
        display: block;
        font-size: 1.25em;
        font-weight: bold;
        margin-bottom: 10px;
      }

      p {
        margin-top: 5px;
        margin-bottom: 5px;
      }
    }
    .content-images {
      width: 40%;
      display: flex;
      flex-direction: column;
      align-items: flex-start;
      justify-content: flex-start;

      padding-left: 10px;

      div {
        margin-bottom: 5px;
      }
    }
  }
`;

const VELogoText = styled.span`
  font-size: 1.5em;
  font-weight: bold;
`;

const ImageRow = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-start;
  a {
    margin-right: 10px;
  }
`;

// Since this component is simple and static, there's no parent container for it.
const AboutTab = () => {
  return (
    <AboutTabContainer>
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
          <span className="section-title">Contributors</span>

          <div className="content-block">
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
              <ul>
                <li>Albert Tian Chen (Vector Engineering, Broad Institute)</li>
                <li>Kevin Altschuler</li>
                <li>
                  Alina Yujia Chan, PhD (Vector Engineering, Broad Institute)
                </li>
                <li>Shing Hei Zhan (University of British Columbia)</li>
                <li>Ben Deverman, PhD (Vector Engineering, Broad Institute)</li>
              </ul>
              <p>
                Contact the authors by email:{' '}
                <a href="mailto:bdeverma@broadinstitute.org">
                  bdeverma@broadinstitute.org
                </a>
                , <a href="mailto:zhan@zoology.ubc.ca">zhan@zoology.ubc.ca</a>,{' '}
                <a href="mailto:alinac@broadinstitute.org">
                  alinac@broadinstitute.org
                </a>
              </p>
              <p>A manuscript for this project is currently being prepared.</p>
            </div>
            <div className="content-images">
              <ImageRow>
                <ExternalLink href="https://vector-engineering.github.io/">
                  <img
                    src="https://storage.googleapis.com/ve-public/covid_ui/assets/img/ve_logo.png"
                    height="40"
                  ></img>
                </ExternalLink>
                <VELogoText>Vector Engineering Lab</VELogoText>
              </ImageRow>
              <ImageRow>
                <ExternalLink href="https://www.broadinstitute.org/">
                  <img
                    src="https://storage.googleapis.com/ve-public/covid_ui/assets/img/BroadLogo_RGB_forDigital.png"
                    height="40"
                    style={{ marginRight: '10px' }}
                  ></img>
                </ExternalLink>
                <ExternalLink href="https://www.broadinstitute.org/stanley">
                  <img
                    src="https://storage.googleapis.com/ve-public/covid_ui/assets/img/StanleyCenterLogo_RGB_forDigital.png"
                    height="50"
                  ></img>
                </ExternalLink>
              </ImageRow>
              <ImageRow>
                <ExternalLink href="https://www.ubc.ca/">
                  <img src={UBCLogo} height="60"></img>
                </ExternalLink>
              </ImageRow>
            </div>
          </div>
        </ContentSection>

        <ContentSection>
          <a id="sequence-data"></a>
          <span className="section-title">Sequence Data</span>

          <div className="content-block">
            <div className="content-text">
              <p>
                We are extremely grateful for{' '}
                <ExternalLink href="https://www.gisaid.org/">
                  GISAID
                </ExternalLink>{' '}
                for sharing all 50K+ SARS-CoV-2 sequences. This project would
                not be possible without their EpiCov™ database, associated
                tooling, and data from originating laboratories. Detailed
                acknowledgements for each dataset can be found on the bottom of
                the page they are displayed in, under the
                &quot;Acknowledgements&quot; table.
              </p>
              <p>
                All data use on COVID CG is subject to the GISAID EpiCov™{' '}
                <ExternalLink href="https://www.gisaid.org/registration/terms-of-use/">
                  Database Access Agreement
                </ExternalLink>
              </p>
            </div>
            <div className="content-images">
              <ExternalLink href="https://www.gisaid.org/">
                <img
                  src="https://storage.googleapis.com/ve-public/covid_ui/assets/img/gisaid.png"
                  height="60"
                ></img>
              </ExternalLink>
            </div>
          </div>
        </ContentSection>

        <ContentSection>
          <a id="genomics-projects"></a>
          <span className="section-title">SARS-CoV-2 Genomics Projects</span>

          <div className="content-block">
            <div className="content-text">
              <p>
                <b>
                  <ExternalLink href="https://nextstrain.org/ncov/global">
                    NextStrain nCoV
                  </ExternalLink>
                  : analysis and visualization of pathogen sequence data
                </b>
              </p>
              <p>
                &quot;Nextstrain is an open-source project to harness the
                scientific and public health potential of pathogen genome data.
                We provide a continually-updated view of publicly available data
                with powerful analytics and visualizations showing pathogen
                evolution and epidemic spread. Our goal is to aid
                epidemiological understanding and improve outbreak
                response.&quot;
              </p>
            </div>
            <div className="content-images">
              <ExternalLink href="https://nextstrain.org/ncov/global">
                <img src={NextstrainLogo} height="70"></img>
              </ExternalLink>
            </div>
          </div>

          <div className="content-block">
            <div className="content-text">
              <p>
                <b>
                  <ExternalLink href="https://covidep.ust.hk/">
                    COVIDep
                  </ExternalLink>
                  : Real-time reporting of vaccine target recommendations for
                  the COVID-19 coronavirus (SARS-CoV-2)
                </b>
              </p>
              <p>
                &quot;COVIDep is developed and maintained by the Signal
                Processing &amp; Computational Biology Lab (
                <ExternalLink href="https://www.mckayspcb.com/">
                  SPCB
                </ExternalLink>
                ) at the Hong Kong University of Science and Technology
                (HKUST).&quot;
              </p>
              <p>
                &quot;COVIDep aims to provide real-time potential vaccine
                targets for SARS-CoV-2. It screens the SARS-derived B cell and T
                cell epitopes (available at{' '}
                <ExternalLink href="https://www.viprbrc.org/brc/home.spg?decorator=corona">
                  VIPR
                </ExternalLink>{' '}
                / <ExternalLink href="http://www.iedb.org/">IEDB</ExternalLink>)
                and identifies those which are highly conserved within the
                available SARS-CoV-2 sequences (continuing to be deposited at{' '}
                <ExternalLink href="https://www.gisaid.org/CoV2020/">
                  GISAID
                </ExternalLink>
                ).&quot;
              </p>
            </div>
            <div className="content-images">
              <ExternalLink href="https://covidep.ust.hk/">
                <img src={HKUSTLogo} height="80"></img>
              </ExternalLink>
            </div>
          </div>

          <div className="content-block">
            <div className="content-text">
              <p>
                <b>
                  <ExternalLink href="https://cov.lanl.gov/content/index">
                    COVID-19 Viral Genome Analysis Pipeline
                  </ExternalLink>
                </b>
              </p>
              <p>
                &quot;This website provides analyses and tools for exploring
                accruing mutations in hCoV-19 (SARS-CoV-2) geographically and
                over time, with an emphasis on the Spike protein, using data
                from GISAID.&quot;
              </p>
              <p>
                &quot;The details of the analyses are described in: Korber et
                al., <i>Cell</i>, June 2020. DOI:{' '}
                <ExternalLink href="https://doi.org/10.1016/j.cell.2020.06.043">
                  10.1016/j.cell.2020.06.043
                </ExternalLink>
                &quot;
              </p>
            </div>
            <div className="content-images">
              <ExternalLink href="https://cov.lanl.gov/content/index">
                <img src={LANLLogo} height="90"></img>
              </ExternalLink>
            </div>
          </div>

          <div className="content-block">
            <div className="content-text">
              <p>
                <b>
                  <ExternalLink href="https://cov-lineages.org/index.html">
                    SARS-CoV-2 lineages
                  </ExternalLink>
                  : A dynamic nomenclature proposal for SARS-CoV-2 lineages to
                  assist genomic epidemiology.
                </b>
              </p>
              <p>
                Rambaut A, Holmes EC, O’Toole Á, Hill V, McCrone JT, Ruis C, du
                Plessis L &amp; Pybus OG (2020) <i>Nature Microbiology</i> DOI:{' '}
                <ExternalLink href="https://doi.org/10.1038/s41564-020-0770-5">
                  10.1038/s41564-020-0770-5
                </ExternalLink>
                .
              </p>
              <p>
                Includes tools such as{' '}
                <ExternalLink href="https://cov-lineages.org/pangolin.html">
                  pangolin
                </ExternalLink>{' '}
                and{' '}
                <ExternalLink href="https://cov-lineages.org/llama.html">
                  llama
                </ExternalLink>
                .
              </p>
            </div>
            <div className="content-images">
              <ExternalLink href="https://cov-lineages.org/index.html">
                <img src={PangolinLogo} height="100"></img>
              </ExternalLink>
            </div>
          </div>

          <div className="content-block">
            <div className="content-text">
              <p>
                <b>
                  <ExternalLink href="https://macman123.shinyapps.io/ugi-scov2-alignment-screen/">
                    SARS-CoV-2 Alignment Screen
                  </ExternalLink>
                </b>
              </p>
              <p>
                &quot;Visualisation tool providing the distribution of SNPs and
                homoplasies across an alignment of 23,090 high coverage,
                complete, SARS-CoV-2 assemblies, downloaded on the 4th of June
                2020. Supporting information is available in{' '}
                <ExternalLink href="https://doi.org/10.1101/2020.05.21.108506">
                  van Dorp &amp; Richard et al. 2020
                </ExternalLink>{' '}
                and{' '}
                <ExternalLink href="https://doi.org/10.1016/j.meegid.2020.104351">
                  van Dorp et al. 2020
                </ExternalLink>
                .&quot;
              </p>
            </div>
            <div className="content-images">
              <ExternalLink href="https://macman123.shinyapps.io/ugi-scov2-alignment-screen/">
                <img src={UCLLogo} width="100"></img>
              </ExternalLink>
            </div>
          </div>

          <div className="content-block">
            <div className="content-text">
              <p>
                <b>
                  <ExternalLink href="https://microreact.org/project/COGconsortium">
                    COG UK Microreact
                  </ExternalLink>
                </b>
              </p>
              <p>
                &quot;From this website are aiming to share the results of
                genomic analysis that feed into weekly updates to the UK
                Government and Public Health Agencies to help guide its
                healthcare strategies in responding to, and minimising, the
                spread of COVID-19 across the UK.&quot;
              </p>
              <p>
                &quot;The Centre for Genomic Pathogen Surveillance maintain a
                Microreact website which permits continuous evaluation of the
                lineages circulating in the UK (currently updated weekly).&quot;
              </p>
            </div>
            <div className="content-images">
              <ExternalLink href="https://microreact.org/project/COGconsortium">
                <img src={COGUKLogo} height="80"></img>
              </ExternalLink>
            </div>
          </div>
        </ContentSection>

        <ContentSection>
          <a id="case-trackers"></a>
          <span className="section-title">COVID-19 Case Trackers</span>

          <div className="content-block">
            <div className="content-text">
              <p>
                <b>
                  <ExternalLink href="https://coronavirus.jhu.edu/">
                    JHU Coronavirus Resource Center
                  </ExternalLink>
                </b>
              </p>
              <p>
                &quot;Johns Hopkins experts in global public health, infectious
                disease, and emergency preparedness have been at the forefront
                of the international response to COVID-19.&quot;
              </p>
              <p>
                &quot;This website is a resource to help advance the
                understanding of the virus, inform the public, and brief
                policymakers in order to guide a response, improve care, and
                save lives.&quot;
              </p>
            </div>
            <div className="content-images">
              <ExternalLink href="https://coronavirus.jhu.edu/">
                <img src={JHULogo} height="60"></img>
              </ExternalLink>
            </div>
          </div>

          <div className="content-block">
            <div className="content-text">
              <p>
                <b>
                  <ExternalLink href="https://91-divoc.com/pages/covid-visualization/">
                    91-DIVOC
                  </ExternalLink>
                  : An interactive visualization of the exponential spread of
                  COVID-19
                </b>
              </p>
              <p>
                &quot;A project to explore the global growth of COVID-19.
                Updated daily.{' '}
                <ExternalLink href="https://91-divoc.com/pages/covid-visualization/overview.html">
                  Overview and motivations
                </ExternalLink>
                .&quot;
              </p>
            </div>
            <div className="content-images"></div>
          </div>

          <div className="content-block">
            <div className="content-text"></div>
            <div className="content-images"></div>
          </div>
        </ContentSection>

        <ContentSection>
          <a id="reference-data"></a>
          <span className="section-title">Reference Sequence</span>

          <div className="content-block">
            <div className="content-text">
              <p>
                All SARS-CoV-2 sequences are aligned to the Wuhan-Hu-1 reference
                sequence (
                <ExternalLink href="https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2">
                  NCBI: NC_045512.2
                </ExternalLink>
                , &quot;WH01&quot; in Lu et al., 2020, DOI:{' '}
                <ExternalLink href="https://doi.org/10.1016/S0140-6736(20)30251-8">
                  10.1016/S0140-6736(20)30251-8
                </ExternalLink>
                ). Gene and protein ORFs were taken from the GenBank page.
              </p>
            </div>
            <div className="content-images"></div>
          </div>
        </ContentSection>
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

        <ContentSection>
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
        </ContentSection>
      </Content>
    </AboutTabContainer>
  );
};

export default AboutTab;
