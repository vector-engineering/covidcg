import React from 'react';
import styled from 'styled-components';
import { observer } from 'mobx-react';
import { FixedSizeList as List } from 'react-window';
import AutoSizer from 'react-virtualized-auto-sizer';
import TableLegendItem from './TableLegendItem';
import { useStores } from '../../stores/connect';

const StyledContainer = styled.div`
  width: 100%;
  height: 100%;
`;

const Columns = styled.div`
  width: 100%;
  display: flex;
  margin: 4px 2px;
`;

const SortArrow = ({ dir }) => {
  if (dir === 'DESC') {
    return <>&#9660;</>;
  } else if (dir === 'ASC') {
    return <>&#9650;</>;
  }
};

const StyledColumnHeader = styled.div`
  cursor: pointer;
  border-bottom: 1px solid #eee;
  width: 50%;
  font-size: 12px;
  padding-left: 2px;
`;

const TableLegend = observer(
  ({
    legendItems,
    updateHoverGroup,
    updateSelectGroup,
    sortColumn,
    sortDir,
    onClickColumnHeader,
  }) => {

    const { configStore } = useStores();

    const Row = ({ index, style }) => {
      const legendItem = legendItems[index];
      return (
        <TableLegendItem
          style={style}
          group={legendItem.group}
          color={legendItem.color}
          updateHoverGroup={updateHoverGroup}
          updateSelectGroup={updateSelectGroup}
          percentage={legendItem.cases_percent}
        />
      );
    };

    const ColumnHeader = ({ columnName, width, children }) => {
      return (
        <StyledColumnHeader onClick={() => onClickColumnHeader({ columnName })}>
          {children} {sortColumn === columnName && <SortArrow dir={sortDir} />}
        </StyledColumnHeader>
      );
    };

    return (
      <StyledContainer>
        <Columns>
          <ColumnHeader columnName="group" width="55%">
            {configStore.getGroupLabel()}
          </ColumnHeader>
          <ColumnHeader columnName="cases_percent" width="45%">
            % Seqs
          </ColumnHeader>
        </Columns>
        <AutoSizer>
          {({ height, width }) => (
            <List
              className="List"
              height={height}
              itemCount={legendItems ? legendItems.length : 0}
              itemSize={35}
              width={width}
            >
              {Row}
            </List>
          )}
        </AutoSizer>
      </StyledContainer>
    );
  }
);

export default TableLegend;
