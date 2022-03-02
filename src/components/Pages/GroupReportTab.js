import React, { useEffect } from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import { config } from '../../config';

import ReactTooltip from 'react-tooltip';

import AccordionWrapper from '../Common/AccordionWrapper';
import ExternalLink from '../Common/ExternalLink';
import GroupTreePlot from '../Vega/GroupTreePlot';
import GroupReportHeader from '../GroupReport/GroupReportHeader';
import MutationList from '../GroupReport/MutationList';
import StructuralViewer from '../GroupReport/GroupStructuralViewer';

import {
  GroupReportTabContainer,
  GroupTreePlotContainer,
  MutationsContainer,
  VOCContainer,
  StructuralViewerContainer,
} from './GroupReportTab.styles';

const GroupReportTab = observer(() => {
  const { groupDataStore } = useStores();

  useEffect(() => {
    ReactTooltip.rebuild();
  }, []);

  const renderHeader = () => {
    return (
      <AccordionWrapper
        title={`${groupDataStore.getActiveGroupTypePrettyName()} Report`}
        defaultCollapsed={false}
        maxHeight={'600px'}
        helpText={
          <ul>
            <li>
              Select {groupDataStore.getActiveGroupTypePrettyName()}s to analyze
              in the plots below as well as the tree on the left.
            </li>
            <li>
              Use the search function to look for specific{' '}
              {groupDataStore.getActiveGroupTypePrettyName()}s
            </li>
            <li>
              Variants of Concern are being most closely monitored, followed by
              Variants of Interest (called Variants Under Investigation by the
              UK HSA).
            </li>
            <li>
              Any other, organization-specific classifications are collectively
              grouped under Other Variants Being Monitored.
            </li>
            <li>
              If a square next to a lineage is colored in, that organization has
              given the lineage whichever classification it is below. Different
              organizations may give the same lineage different classifications.
            </li>
          </ul>
        }
      >
        {config.virus === 'sars2' && <GroupReportHeader />}
      </AccordionWrapper>
    );
  };

  const renderMutationList = () => {
    return <MutationList />;
  };

  const renderStructuralViewer = () => {
    return (
      <AccordionWrapper
        title={`Structural Viewer`}
        defaultCollapsed={false}
        maxHeight={'1000px'}
        helpText={
          <ul>
            <li>
              Molecule visualizations are provided by{' '}
              <ExternalLink href="https://www.litemol.org/">
                LiteMol
              </ExternalLink>
              . Click the &quot;?&quot; button for control help.
            </li>
            <li>
              Structures are downloaded from the{' '}
              <ExternalLink href="https://www.rcsb.org/">RCSB PDB</ExternalLink>
            </li>
            <li>
              Residues are colored by the mutation frequencies in the selected
              Protein, for the selected{' '}
              {groupDataStore.getGroupMutationTypePrettyName()}, projected onto
              the given PDB ID
            </li>
            <li>
              It is up to the user to ensure that the given Protein and its
              mutations match up with the provided PDB ID. We do not check for
              compatibility, so it is possible to, for example, erroneously map
              nsp12 mutations onto a Spike structure.
            </li>
            <li>
              More structures available at the{' '}
              <ExternalLink href="https://www.rcsb.org/news?year=2020&article=5e74d55d2d410731e9944f52&feature=true">
                RCSB COVID-19 resource page
              </ExternalLink>
            </li>
          </ul>
        }
      >
        <StructuralViewer />
      </AccordionWrapper>
    );
  };

  return (
    <GroupReportTabContainer>
      <GroupTreePlotContainer>
        {config.virus === 'sars2' && <GroupTreePlot width={250} />}
      </GroupTreePlotContainer>
      <MutationsContainer>{renderMutationList()}</MutationsContainer>
      <VOCContainer>{renderHeader()}</VOCContainer>
      <StructuralViewerContainer>
        {renderStructuralViewer()}
      </StructuralViewerContainer>
    </GroupReportTabContainer>
  );
});

export default GroupReportTab;
