import React from 'react';
import styled from 'styled-components';

import ExternalLink from '../Common/ExternalLink';
import VELogo from '../../assets/images/VE_logo_new.png';

import { config } from '../../config';

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
        <ContentSection>
          <a id="about"></a>
          <span className="section-title">About COVID CG</span>

          <div className="content-text">
            <p>
              In the first year since COVID CG was launched, the site has been
              visited more than 47,900 times by users from 180 countries. We
              have fielded inquiries from AP News, CNN, WSJ, and scientists from
              the European Commission, Canadian government, and Australian
              government who are using our site to track genomic variants and
              global sequencing efforts. Thanks to the GISAID database, the
              covidcg.org site is helping to make an impact on the development
              of Covid-19 diagnostics and therapeutics in different countries.
            </p>
            <p>
              Our new site, with significant behind the scenes enhancements,
              allows users to quickly analyze the 2.2 million and growing
              SARS-CoV-2 genomes deposited in GISAID. We are also excited to
              present our newly released Lineage Reports page.
            </p>
            <p>
              <b>
                On the COVID CG’s Lineage Reports page, you can easily learn
                more about Variants of Concern (VOC) and Variants of Interest
                (VOI) using these new features:
              </b>
              <ul>
                <li>
                  <b>Phylogenetic Time-Scaled Tree:</b> How are the different
                  SARS-CoV-2 lineages and variants related to each other? This
                  feature was developed based on open source code shared by Art
                  Poon who also runs a GISAID-powered{' '}
                  <ExternalLink href="https://filogeneti.ca/covizu/">
                    CoVizu
                  </ExternalLink>
                  resource.
                </li>
                <li>
                  <b>Mutation Type Heat Map:</b> What are the mutations that
                  distinguish my lineage or variant from the original SARS-CoV-2
                  genome? How frequently do specific mutations occur in that
                  lineage or variant as compared to other lineages? This
                  analysis covers all SARS-CoV-2 genes, not just the spike.
                </li>
                <li>
                  <b>Structural Viewer:</b> For Spike mutations in my lineage or
                  variant of interest, how do these map onto the protein
                  structure of the SARS-CoV-2 Spike?
                </li>
              </ul>
            </p>
            <p>
              <b>
                Other interactive features have been upgraded so that users can
                rapidly visualize:
              </b>
              <ul>
                <li>
                  The emergence of new lineages and variants in a location and
                  date-specific manner.
                </li>
                <li>
                  Genome-wide entropy and mutations that co-occur with your
                  mutation of interest.
                </li>
              </ul>
            </p>
            <p>
              COVID CG continues to be a free, public resource, supported by{' '}
              <ExternalLink href="https://giving.broadinstitute.org/broadignite">
                Broad Ignite
              </ExternalLink>
              and a research collaboration with AstraZeneca.
            </p>
          </div>
        </ContentSection>
        <ContentSection>
          <a id="contributors"></a>
          <span className="section-title">COVID-CG is developed by</span>

          <div className="content-text">
            <p>
              COVID CoV Genomics (CG) was developed in the
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
              <b>David Favela</b> is an associate software engineer in the
              Deverman group.
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
              <ExternalLink href="https://www.broadinstitute.org/bios/ben-deverman">
                Dr. Benjamin Deverman
              </ExternalLink>{' '}
              is the director of the vector engineering research group at the
              Stanley Center for Psychiatric Research at the Broad Institute of
              MIT and Harvard, where he is also an Institute Scientist. The
              Deverman Lab develops innovative gene delivery solutions to expand
              the impact of gene therapy and for studying the central nervous
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

        {config['show_logos']['GISAID'] && (
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
        )}

        <ContentSection>
          <a id="citing-covid-cg"></a>
          <span className="section-title">Citing COVID-19 CG:</span>
          <div className="content-text">
            <p>COVID-19 CG is published in eLife:</p>
            <p>
              Chen AT, Altschuler K, Chan AY, Zhan SH, Deverman BE (2021).
              COVID-19 CG enables SARS-CoV-2 mutation and lineage tracking by
              locations and dates of interest. <i>eLife</i>. DOI:{' '}
              <ExternalLink href="https://doi.org/10.7554/eLife.63409">
                10.7554/eLife.63409
              </ExternalLink>
            </p>
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
      </Content>
    </TabContainer>
  );
};

export default AboutTab;
