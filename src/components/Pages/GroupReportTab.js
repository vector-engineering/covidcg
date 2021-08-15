import React, { useEffect } from 'react';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

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
              Notable lineages, as defined by the{' '}
              <ExternalLink href="https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-info.html">
                CDC
              </ExternalLink>
              , can be selected from the displayed list
            </li>
          </ul>
        }
      >
        <GroupReportHeader />
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
              {groupDataStore.getGroupSnvTypePrettyName()}, projected onto the
              given PDB ID
            </li>
            <li>
              It is up to the user to ensure that the given Protein and its SNVs
              match up with the provided PDB ID. We do not check for
              compatibility, so it is possible to, for example, erroneously map
              nsp12 SNVs onto a Spike structure.
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
        <GroupTreePlot width={250} />
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
