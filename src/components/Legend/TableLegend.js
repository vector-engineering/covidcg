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
  width: ${({ width }) => width};
  font-size: 12px;
  padding-left: 2px;
`;

/*

*/

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
          percentage={legendItem.percent}
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
        <Columns>
          <ColumnHeader columnName="group" width="55%">
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
          <ColumnHeader columnName="percent" width="45%">
            % Seqs
          </ColumnHeader>
        </Columns>
        <AutoSizer>
          {({ height, width }) => (
            <List
              className="List"
              height={height - 30}
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
