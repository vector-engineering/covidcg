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
  border-left: 1px solid black;
  border-right: 1px solid black;
  cursor: col-resize;
  position: fixed;
  right: 400px;
  top: 0;
  z-index: 2001;
  box-shadow: 0px 0px 63px -13px rgba(0, 0, 0, 1);
  display: flex;
  align-items: center;
  justify-content: center;
  font-size: 24px;
  color: black;
  background-color: white;
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
    width: `${props.width - 40}px`,
  },
}))`
  background-color: white;
  height: 100%;
  padding-right: 20px;
  padding-left: 20px;
  overflow-x: hidden;
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
  flex-direction: column;
  width: 100%;
  height: 100%;
  overflow-y: scroll;
`;

const CloseOpenButton = styled.button.attrs((props) => ({
  style: {
    right: `${props.right}px`,
  },
}))`
  position: absolute;
  top: 5px;
  border: 0px;
  background-color: #eee;
  font-size: 32px;
  cursor: pointer;
  &:hover {
    background-color: #bbb;
  }
`;

const LiteMolBlock = styled.div`
  width: 100%;
  height: 100%;
`;

const SideBar = observer(() => {
  const [widthDelta, setWidthDelta] = useState(400);
  const { uiStore } = useStores();

  const onOpenSidebar = () => {
    setWidthDelta(400);
    uiStore.setSidebarOpen();
  };
  const onCloseSidebar = () => {
    setWidthDelta(400);
    uiStore.setSidebarClosed();
  };

  if (!uiStore.sidebarOpen) {
    return (
      <SideBarClosedContainer>
        <CloseOpenButton right={0} onClick={onOpenSidebar}>
          {'<'}
        </CloseOpenButton>
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
        <ClickDragBorder>⫶</ClickDragBorder>
      </Draggable>
      <SidebarContainer>
        <SidebarContent width={widthDelta}>
          <CloseOpenButton right={widthDelta + 20} onClick={onCloseSidebar}>
            &gt;
          </CloseOpenButton>
          <ViewersContainer>
            {uiStore.sidebarSelectedGroupKeys.map((groupKey) => (
              <LiteMolBlock key={groupKey}>
                <LiteMolBlockHeader>
                  {groupKey}
                  <button
                    onClick={() => uiStore.onRemoveGroupFromSidebar(groupKey)}
                  >
                    x
                  </button>
                </LiteMolBlockHeader>
                <LiteMolViewer key={Math.random()} />
              </LiteMolBlock>
            ))}
          </ViewersContainer>
        </SidebarContent>
      </SidebarContainer>
    </>
  );
});

export default SideBar;
