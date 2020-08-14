import React from 'react';

import ExternalLink from '../Common/ExternalLink';

import {
  TabContainer,
  Content,
  ContentSection,
  // ImageRow,
} from './TextTab.styles';

import WorkflowImage from '../../assets/images/Fig2_workflow_V4.svg';

const MethodologyTab = () => {
  return (
    <TabContainer>
      <Content style={{ width: '100%' }}>
        <ContentSection>
          <a id="workflow-fig"></a>
          <span className="section-title">Workflow</span>
          <div style={{ height: '10px' }} />
          <img style={{ width: '100%' }} src={WorkflowImage} />
        </ContentSection>

        <ContentSection>
          <a id="preprocessing"></a>
          <span className="section-title">Sequence Preprocessing</span>

          <div className="content-block">
            <div className="content-text">
              <p>
                Sequences meeting any of the following criteria are filtered
                out:
              </p>
              <ul>
                <li>
                  Present on the{' '}
                  <ExternalLink href="https://github.com/nextstrain/ncov/blob/master/defaults/exclude.txt">
                    Nextstrain exclusion list
                  </ExternalLink>
                </li>
                <li>
                  Isolates from non-humans (animals, enviornmental samples, etc)
                </li>
                <li>Less than 29700 bases</li>
                <li>&gt; 5% ambiguous base calls (N)</li>
              </ul>
            </div>
            <div className="content-images"></div>
          </div>
        </ContentSection>

        <ContentSection>
          <a id="metadata-cleaning"></a>
          <span className="section-title">Metadata Cleaning</span>

          <div className="content-block">
            <div className="content-text">
              <p>
                Metadata cleaning aims to preserve the original intent of the
                authors and data submitters while presenting simpler and unified
                versions to end users.
              </p>
              <p>
                Sequence metadata is cleaned to remove obvious typos, and to
                unify labels with the same meaning, e.g., &quot;MinION&quot; and
                &quot;Nanopore MinION&quot;.
              </p>
              {/* <p>
                Patient ages are transformed into ranges, in order to work
                across the many levels of detail present in patient age data.
                For example &quot;65&quot; will be interpreted as [65, 66)
              </p> */}
              <p>
                Location metadata is cleaned with the goal of simplifying the
                location selector in the sidebar. Locations with excessive
                children are collapsed to the nearest upper hierarchical
                grouping. E.g., if a state has individual data for 200+ towns,
                these towns will be collapsed to the county level in order to
                facilitate easier data browsing. Typos and clear similarities
                are also unified to prevent duplicate locations
              </p>
            </div>
          </div>
        </ContentSection>

        <ContentSection>
          <a id="alignment"></a>
          <span className="section-title">SNV Assignments</span>

          <div className="content-block">
            <div className="content-text">
              <p>
                SNVs at the nucleotide and amino acid level will be determined
                by aligning each sequence to the{' '}
                <ExternalLink href="https://www.ncbi.nlm.nih.gov/nuccore/MN996528">
                  WIV04 reference sequence
                </ExternalLink>{' '}
                (WIV04 is 100% identical to a high quality December, 2019
                isolate, Wuhan-Hu-1/
                <ExternalLink href="https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2">
                  NC_045512.2
                </ExternalLink>{' '}
                used by NextStrain) using bowtie2. Spurious SNVs and probable
                sequencing errors are filtered out prior to downstream analysis.
                SNVs involving ambiguous base calls are ignored.
              </p>
            </div>
          </div>
        </ContentSection>

        <ContentSection>
          <a id="lineage-analysis"></a>
          <span className="section-title">Lineage/Clade Analysis</span>

          <div className="content-block">
            <div className="content-text">
              <p>
                Viral lineages (as defined by the{' '}
                <ExternalLink href="https://github.com/cov-lineages/pangolin">
                  pangolin
                </ExternalLink>{' '}
                tool), and clades will be provided by GISAID. In accordance with
                pangolin, SNVs present in &gt;90% of sequences within each
                lineage/clade will be assigned as lineage/clade-defining SNVs.
              </p>
            </div>
          </div>
        </ContentSection>
      </Content>
    </TabContainer>
  );
};

export default MethodologyTab;
