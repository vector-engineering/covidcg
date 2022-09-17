import React from 'react';

import { config } from '../../config';

import ExternalLink from '../Common/ExternalLink';

import {
  AcknowledgementFooterContainer,
  SectionTitle,
  ContentText,
} from './AcknowledgementFooter.styles';

const AcknowledgementFooter = ({ ...props }) => {
  return (
    <AcknowledgementFooterContainer {...props}>
      {(config['data_provider'] === 'GISAID' ||
        config['postgres_db'] == 'cg_genbank_dev') && (
        <>
          <SectionTitle>Data enabling {config.site_title}</SectionTitle>
          <ContentText>
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
          </ContentText>
        </>
      )}

      <SectionTitle>Citing {config.site_title}:</SectionTitle>
      <ContentText>
        <p>COVID-19 CG is published in eLife:</p>
        <p>
          Chen AT, Altschuler K, Chan AY, Zhan SH, Deverman BE (2021). COVID-19
          CG enables SARS-CoV-2 mutation and lineage tracking by locations and
          dates of interest. <i>eLife</i>. DOI:{' '}
          <ExternalLink href="https://doi.org/10.7554/eLife.63409">
            10.7554/eLife.63409
          </ExternalLink>
        </p>
        <p>
          Users are encouraged to share, download, and further analyze data from
          this site. Plots can be downloaded as PNG or SVG files, and the data
          powering the plots and tables can be downloaded as well. Please
          attribute any data/images to{' '}
          <ExternalLink href={config.prod_hostname} title={config.site_title}>
            {config.prod_hostname}
          </ExternalLink>
          .
        </p>
        {(config['data_provider'] === 'GISAID' ||
          config['postgres_db'] == 'cg_genbank_dev') && (
          <>
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
            <p>
              All data in GISAID are subject to GISAID's{' '}
              <ExternalLink href="https://www.gisaid.org/registration/terms-of-use/">
                Terms and Conditions
              </ExternalLink>
              .
            </p>
          </>
        )}
      </ContentText>
    </AcknowledgementFooterContainer>
  );
};

export default AcknowledgementFooter;
