import React, { useEffect, useState, useRef } from 'react';
import styled from 'styled-components';

const CellContainer = styled.div`
  display: flex;
  flex-direction: row;
`;

const LiteMolViewer = () => {
  const [isOpen, setIsOpen] = useState(false);

  const uniqueId = useRef(`litemolcell-${new Date().getTime()}`);

  useEffect(() => {
    const litemoldiv = document.querySelector(`#${uniqueId.current}`);
    if (isOpen) {
      window.angular.element(litemoldiv).ready(function () {
        window.angular.bootstrap(litemoldiv, ['pdb.litemol']);
      });
    }
  }, [isOpen]);

  return (
    <CellContainer>
      {isOpen && (
        <div style={{ position: 'relative', height: '200px', width: 400 }}>
          <div
            source-format="sdf"
            source-url="http://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/files/sdf/ATP.sdf"
            id={uniqueId.current}
            className="pdb-lite-mol"
            pdb-id="'1cbs'"
          ></div>
        </div>
      )}
      <button
        onClick={() => {
          setIsOpen(!isOpen);
        }}
      >
        open/close
      </button>
    </CellContainer>
  );
};

export default LiteMolViewer;
