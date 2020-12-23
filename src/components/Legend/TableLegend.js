import React from 'react';
import { observer } from 'mobx-react';
import _ from 'underscore';
import { FixedSizeList as List } from 'react-window';
import AutoSizer from 'react-virtualized-auto-sizer';
import TableLegendItem from './TableLegendItem';

const TableLegend = observer(
  ({ legendItems, updateHoverGroup, updateSelectGroup }) => {
    const Row = ({ index, style }) => {
      const legendItem = legendItems[index];
      return (
        <TableLegendItem
          style={style}
          group={legendItem.group}
          color={legendItem.color}
          updateHoverGroup={updateHoverGroup}
          updateSelectGroup={updateSelectGroup}
        />
      );
    };

    return (
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
    );
  }
);

export default TableLegend;
