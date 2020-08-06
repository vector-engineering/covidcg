import React from 'react';

import ExternalLink from '../Common/ExternalLink';

import NextstrainLogo from '../../assets/images/nextstrain_logo.png';
import HKUSTLogo from '../../assets/images/HKUST-original_0.svg';
import LANLLogo from '../../assets/images/lanl_logo.svg';
import PangolinLogo from '../../assets/images/pangolin_logo.png';
import UCLLogo from '../../assets/images/university-college-london-ucl-vector-logo.svg';
import COGUKLogo from '../../assets/images/logo-cog-uk.png';
import JHULogo from '../../assets/images/jhu_logo.jpg';

import {
  TabContainer,
  Content,
  ContentSection,
  ImageRow,
} from './TextTab.styles';

const RelatedProjectsTab = () => {
  return (
    <TabContainer>
      <Content>
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
      </Content>
    </TabContainer>
  );
};

export default RelatedProjectsTab;
