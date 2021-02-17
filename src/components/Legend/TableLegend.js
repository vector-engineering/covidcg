import React from 'react';
import styled from 'styled-components';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import AutoSizer from 'react-virtualized-auto-sizer';
import { FixedSizeList as List } from 'react-window';
import QuestionButton from '../Buttons/QuestionButton';
import ReactTooltip from 'react-tooltip';
import TableLegendItem from './TableLegendItem';

const StyledContainer = styled.div`
  width: 100%;
  height: 100%;
`;

const Header = styled.div`
  height: 48px;
  border-bottom: 1px solid #ccc;
`;

const HeaderRow = styled.div`
  width: 100%;
  height: 50%;
  display: flex;
  flex-direction: row;
  align-items: center;
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
  width: ${({ width }) => width};
  font-size: 12px;
  padding: 0px 3px;
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

    const maxCounts = legendItems.reduce((prev, cur) => {
      return cur.counts > prev ? cur.counts : prev;
    }, 0);

    const Row = ({ index, style }) => {
      const legendItem = legendItems[index];
      return (
        <TableLegendItem
          style={style}
          item={legendItem}
          maxCounts={maxCounts}
          updateHoverGroup={updateHoverGroup}
          updateSelectGroup={updateSelectGroup}
        />
      );
    };
    Row.propTypes = {
      index: PropTypes.number,
      style: PropTypes.object,
    };

    const ColumnHeader = ({ columnName, width, children }) => {
      return (
        <StyledColumnHeader
          onClick={() => onClickColumnHeader({ columnName })}
          width={width}
        >
          {children} {sortColumn === columnName && <SortArrow dir={sortDir} />}
        </StyledColumnHeader>
      );
    };
    ColumnHeader.propTypes = {
      columnName: PropTypes.string,
      width: PropTypes.string,
      children: PropTypes.oneOfType([
        PropTypes.arrayOf(PropTypes.node),
        PropTypes.node,
      ]),
    };
    ColumnHeader.defaultProps = {
      width: '50%',
    };

    return (
      <StyledContainer>
        <Header>
          <HeaderRow>
            <ColumnHeader columnName="group" width="100%">
              <ReactTooltip
                className="legend-sidebar-tooltip"
                id="legend-sidebar-tooltip"
                type="light"
                effect="solid"
                border={true}
                borderColor="#888"
              />
              {configStore.getGroupLabel()}
              <QuestionButton
                data-tip={`
              <ul>
                <li>Items in the legend represent <b>${configStore.getGroupLabel()}s</b>.
                </li>
                <li>
                  Click to select one, or hold Shift and click to select
                  multiple ${configStore.getGroupLabel()}s.
                </li>
                <li>
                  Selected ${configStore.getGroupLabel()}s will be highlighted in the
                  plots and table below.
                </li>
              </ul>
              `}
                data-html="true"
                data-for="legend-sidebar-tooltip"
              />
            </ColumnHeader>
          </HeaderRow>
          <HeaderRow>
            <ColumnHeader columnName="spacer" width="40%"></ColumnHeader>
            <ColumnHeader columnName="counts" width="30%">
              #
            </ColumnHeader>
            <ColumnHeader columnName="percent" width="30%">
              %
            </ColumnHeader>
          </HeaderRow>
        </Header>
        <AutoSizer>
          {({ height, width }) => (
            <List
              className="List"
              height={height - 48}
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
