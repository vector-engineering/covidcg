import React, { useEffect, useRef } from 'react';

const LiteMolViewer = () => {
  const uniqueId = useRef(`litemolcell-${new Date().getTime()}`);

  useEffect(() => {
    const litemoldiv = document.querySelector(`#${uniqueId.current}`);
    window.angular.element(litemoldiv).ready(function () {
      window.angular.bootstrap(litemoldiv, ['pdb.litemol']);
    });
  });

  return (
    <div style={{ position: 'relative', height: '200px', width: 400 }}>
      <div
        source-format="sdf"
        source-url="http://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/files/sdf/ATP.sdf"
        id={uniqueId.current}
        className="pdb-lite-mol"
        pdb-id="'1cbs'"
        hide-controls="true"
      ></div>
    </div>
  );
};

export default LiteMolViewer;
