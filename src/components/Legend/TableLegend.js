import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import AutoSizer from 'react-virtualized-auto-sizer';
import { FixedSizeList as List } from 'react-window';
import ReactTooltip from 'react-tooltip';

import QuestionButton from '../Buttons/QuestionButton';
import TableLegendItem from './TableLegendItem';

import { SORT_DIRECTIONS } from '../../constants/defs.json';
import { COLUMN_NAMES } from './legendUtils';

import {
  StyledContainer,
  Header,
  HeaderRow,
  StyledColumnHeader,
} from './TableLegend.styles';

const SortArrow = ({ dir }) => {
  if (dir === SORT_DIRECTIONS.SORT_DESC) {
    return <>&#9660;</>;
  } else if (dir === SORT_DIRECTIONS.SORT_ASC) {
    return <>&#9650;</>;
  }
};

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
            <ColumnHeader columnName={COLUMN_NAMES.GROUP} width="100%">
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
            <ColumnHeader columnName={COLUMN_NAMES.COUNTS} width="30%">
              #
            </ColumnHeader>
            <ColumnHeader columnName={COLUMN_NAMES.PERCENT} width="30%">
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
