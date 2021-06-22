import React from 'react';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';
import useDimensions from 'react-use-dimensions';

import AccordionWrapper from '../Common/AccordionWrapper';
import ExternalLink from '../Common/ExternalLink';
import GroupTreePlot from '../Vega/GroupTreePlot';
import GroupReportHeader from '../GroupReport/GroupReportHeader';
import MutationList from '../GroupReport/MutationList';
import StructuralViewer from '../LiteMol/StructuralViewer';

const GroupReportTabContainer = styled.div`
  display: flex;
  flex-direction: row;
`;
const GroupTreePlotContainer = styled.div``;

const MainContainer = styled.div`
  flex-grow: 1;
`;

const GroupReportTab = observer(() => {
  const { groupDataStore } = useStores();
  const [ref, { width }] = useDimensions();

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
    return (
      <AccordionWrapper
        title={`${groupDataStore.getGroupSnvTypePrettyName()} per ${groupDataStore.getActiveGroupTypePrettyName()}`}
        defaultCollapsed={false}
        maxHeight={'600px'}
        helpText={
          <ul>
            <li>
              Note: We define ORF1a and ORF1ab as separate genes. In
              &quot;NT&quot; or &quot;Gene AA&quot; mode, an SNV in ORF1a will
              also be listed in ORF1ab. Switch to &quot;Protein AA&quot; mode to
              see SNVs in the context of proteins (i.e., NSPs)
            </li>
          </ul>
        }
      >
        <MutationList />
      </AccordionWrapper>
    );
  };

  const renderStructuralViewer = () => {
    return (
      <AccordionWrapper
        title={`Structural Viewer`}
        defaultCollapsed={false}
        maxHeight={'1000px'}
        helpText={<ul></ul>}
      >
        <StructuralViewer />
      </AccordionWrapper>
    );
  };

  return (
    <GroupReportTabContainer ref={ref}>
      <GroupTreePlotContainer>
        <GroupTreePlot width={300} />
      </GroupTreePlotContainer>
      <MainContainer>
        {renderHeader()}
        {renderMutationList()}
        {renderStructuralViewer()}
      </MainContainer>
    </GroupReportTabContainer>
  );
});

export default GroupReportTab;
