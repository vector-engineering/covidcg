import React, { useEffect, useRef } from 'react';

const LiteMolViewer = () => {
  const uniqueId = useRef(`litemolcell-${new Date().getTime()}`);
  const liteMolScope = useRef();

  const bindPdbComponentScope = function (element) {
    return window.angular.element(element).isolateScope();
  };

  const highlight = () => {
    //hightlighting

    let selectionDetails = {
      entity_id: '1',
      struct_asym_id: 'A',
      start_residue_number: 100,
      end_residue_number: 300,
    };
    liteMolScope.current.LiteMolComponent.highlightOn(selectionDetails);
  };

  useEffect(() => {
    const litemoldiv = document.querySelector(`#${uniqueId.current}`);

    window.angular.element(litemoldiv).ready(function () {
      window.angular.bootstrap(litemoldiv, ['pdb.litemol']);
      liteMolScope.current = bindPdbComponentScope(litemoldiv);
    });
  });

  return (
    <div
      style={{
        position: 'relative',
        height: '100%',
        width: '100%',
      }}
    >
      <button onClick={highlight}>button</button>
      <div
        style={{
          position: 'relative',
          height: '100%',
          width: '100%',
        }}
      >
        <div
          source-format="pdb"
          source-url="https://files.rcsb.org/download/6X2A.pdb"
          id={uniqueId.current}
          className="pdb-lite-mol"
          pdb-id="'6x2a'"
          hide-controls="true"
        ></div>
      </div>
    </div>
  );
};

export default LiteMolViewer;
