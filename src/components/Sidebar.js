import React, { useState } from 'react';
import { observer } from 'mobx-react';

import styled from 'styled-components';
import Draggable from 'react-draggable'; // <DraggableCore>
import { useStores } from '../stores/connect';
import LiteMolViewer from './LiteMolViewer';

function clamp(num, min, max) {
  return num <= min ? min : num >= max ? max : num;
}

const ClickDragBorder = styled.div`
  width: 10px;
  height: 100%;
  background-color: gray;
  border-left: 2px solid black;
  border-right: 2px solid black;
  cursor: col-resize;
  position: fixed;
  right: 400px;
  top: 0;
  z-index: 2001;
  box-shadow: 0px 0px 63px -13px rgba(0, 0, 0, 1);
`;

const SidebarContainer = styled.div`
  position: fixed;
  right: 0;
  top: 0;
  height: 100vh;
  z-index: 2000;
`;

const SidebarContent = styled.div.attrs((props) => ({
  style: {
    width: `${props.width}px`,
  },
}))`
  background-color: white;
  height: 100%;
`;

const SideBarClosedContainer = styled.div`
  position: fixed;
  top: 200;
  right: 0;
`;

const LiteMolBlockHeader = styled.div`
  display: flex;
  flex-direction: row;
`;

const ViewersContainer = styled.div`
  display: flex;
  width: 100%;
  flex-wrap: wrap;
`;

const SideBar = observer(() => {
  const [widthDelta, setWidthDelta] = useState(400);
  const { uiStore } = useStores();

  if (!uiStore.sidebarOpen) {
    return (
      <SideBarClosedContainer>
        <button onClick={uiStore.openSidebar}>open</button>
        {uiStore.sidebarSelectedGroupKeys.map((groupKey) => (
          <div key={groupKey}>{groupKey}</div>
        ))}
      </SideBarClosedContainer>
    );
  }

  return (
    <>
      <Draggable
        axis="x"
        defaultClassName="DragHandle"
        defaultClassNameDragging="DragHandleActive"
        bounds={{ top: 0, left: -600, right: 0, bottom: 0 }}
        onDrag={(e, data) => {
          setWidthDelta(clamp(widthDelta - data.deltaX, 400, 1000));
        }}
      >
        <ClickDragBorder />
      </Draggable>
      <SidebarContainer>
        <SidebarContent width={widthDelta}>
          <button onClick={uiStore.closeSidebar}>close</button>
          <ViewersContainer>
            {uiStore.sidebarSelectedGroupKeys.map((groupKey) => (
              <div key={groupKey}>
                <LiteMolBlockHeader>
                  {groupKey}
                  <button
                    onClick={() => uiStore.onRemoveGroupFromSidebar(groupKey)}
                  >
                    x
                  </button>
                </LiteMolBlockHeader>
                <LiteMolViewer />
              </div>
            ))}
          </ViewersContainer>
        </SidebarContent>
      </SidebarContainer>
    </>
  );
});

export default SideBar;
