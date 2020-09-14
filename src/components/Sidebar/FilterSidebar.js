import React from 'react';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';

import ReactTooltip from 'react-tooltip';
import QuestionButton from '../Buttons/QuestionButton';
import Header from './Header';
import GroupBySelect from './GroupBySelect';
import CoordinateSelect from './CoordinateSelect';
import MetaFieldSelect from './MetaFieldSelect';
import DropdownContainer from './DropdownContainer';
import FilterDataIntoOther from './FilterDataIntoOther';
import SidebarAccordionWrapper from './SidebarAccordionWrapper';

const FilterSidebarContainer = styled.div`
  position: fixed;
  top: 0;
  width: 299px;

  background-color: #f8f8f8;
  //padding-right: 10px;
  padding-bottom: 15px;
  border-right: 1px solid #aaa;
  display: flex;
  flex-direction: column;
  height: 100vh;
  overflow-y: hidden;

  .filter-sidebar-tooltip {
    background-color: #fff;
    font-weight: normal;
    p {
      margin-top: 2px;
      margin-bottom: 2px;
    }
  }
`;

const FilterSidebar = observer(() => {
  const { configStore } = useStores();

  return (
    <FilterSidebarContainer>
      <ReactTooltip
        className="filter-sidebar-tooltip"
        id="tooltip-filter-sidebar"
        type="light"
        effect="solid"
        border={true}
        borderColor="#888"
      />
      <Header />
      <GroupBySelect />
      <SidebarAccordionWrapper
        title={
          <div>
            Collapse low frequency data
            <QuestionButton
              data-tip={`<p>${configStore.getGroupLabel()}s that do not meet the following criteria will be grouped into the "Other" group. This is done to increase performance in the app</p><p>Including more groups gives more detail into the data, but may come at the cost of app performance.</p>`}
              data-html={true}
              data-for="tooltip-filter-sidebar"
            />
          </div>
        }
        defaultCollapsed={true}
        maxHeight={'250px'}
      >
        <FilterDataIntoOther />
      </SidebarAccordionWrapper>
      <SidebarAccordionWrapper
        title="Genomic coordinates"
        defaultCollapsed={false}
        maxHeight={'420px'}
      >
        <CoordinateSelect />
      </SidebarAccordionWrapper>
      <SidebarAccordionWrapper
        title={
          <div>
            Filter by metadata (advanced)
            <QuestionButton
              data-tip='<p>By default, no filtering is applied on sequence metadata (Default is select all)</p><p>Metadata is dependent on the data submitter, so many fields may be missing and marked as "Unknown".</p><p>Metadata filters are shown in the format "[a &gt; b]", <br/>where "a" is the initial number of sequences matching that metadata field,<br/>and "b" is the number of sequences matching that metadata field after all current metadata filters are applied.</p>'
              data-html={true}
              data-for="tooltip-filter-sidebar"
            />
          </div>
        }
        defaultCollapsed={true}
        maxHeight={'240px'}
        allowOverflow={true}
      >
        <MetaFieldSelect />
      </SidebarAccordionWrapper>

      {/*<SidebarAccordionWrapper
      title="Selected locations"
      defaultCollapsed={false}
    >
      <DropdownContainer />
    </SidebarAccordionWrapper>*/}
      <DropdownContainer />
    </FilterSidebarContainer>
  );
});

export default FilterSidebar;
