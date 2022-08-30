import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import AutoSizer from 'react-virtualized-auto-sizer';
import { FixedSizeList as List } from 'react-window';

import ReactTooltip from 'react-tooltip';
import QuestionButton from '../Buttons/QuestionButton';
import TableLegendItem from './TableLegendItem';

import { SORT_DIRECTIONS, GROUP_MUTATION } from '../../constants/defs.json';
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
    const { configStore, plotSettingsStore } = useStores();

    const toggleLegendAdjustPartialSequences = (event) => {
      plotSettingsStore.setLegendAdjustPartialSequences(event.target.checked);
    };

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

    // Header height offset, in pixels
    // Each header row is 24 px tall
    let headerHeightOffset = 48;
    // If we're in mutation mode, show the partial sequences toggle
    const showPartialSequencesToggle = configStore.groupKey === GROUP_MUTATION;
    if (showPartialSequencesToggle) {
      headerHeightOffset += 24;
    }

    return (
      <StyledContainer>
        <Header>
          {showPartialSequencesToggle && (
            <HeaderRow>
              <input
                id="legend-adjust-partial-sequences-checkbox"
                type="checkbox"
                checked={plotSettingsStore.legendAdjustPartialSequences}
                onChange={toggleLegendAdjustPartialSequences}
              ></input>
              <label htmlFor="legend-adjust-partial-sequences-checkbox">
                Adjust for coverage
              </label>
              <QuestionButton
                data-tip={`
              <ul>
                <li>
                  Some isolates are not sequenced across the entire genome, 
                  i.e., partial sequences of just one gene
                </li>
                <li>
                  This option adjusts the percentages to reflect the coverage 
                  of the genome at each position
                </li>
              </ul>
              `}
                data-html="true"
                data-for="legend-sidebar-tooltip"
              />
            </HeaderRow>
          )}
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
              height={height - headerHeightOffset} // 24 px per header row * 3 rows
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
