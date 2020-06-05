import React, { useEffect, useState } from 'react';
import styled from 'styled-components';

const CellContainer = styled.div`
  display: 'flex';
  flex-direction: 'row';
`;

const LiteMolCell = () => {
  const [isOpen, setIsOpen] = useState(false);

  useEffect(() => {
    const litemoldiv = document.querySelector('#litemol');
    if (isOpen) {
      window.angular.element(litemoldiv).ready(function () {
        console.log(window.angular.bootstrap(litemoldiv, ['pdb.litemol']));
      });
    }
  }, [isOpen]);

  return (
    <CellContainer>
      {isOpen && (
        <div>
          yo
          <div style={{ position: 'relative', height: '200px', width: 400 }}>
            <div
              source-format="sdf"
              source-url="http://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/files/sdf/ATP.sdf"
              id="litemol"
              className="pdb-lite-mol"
              pdb-id="'1cbs'"
            ></div>
          </div>
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

export default LiteMolCell;
