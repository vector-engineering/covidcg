import React, { useState, useEffect, useRef } from 'react';
import './../../styles/litemol.min.css';

import LiteMol from 'litemol';

const Rx = LiteMol.Core.Rx;
const Transformer = LiteMol.Bootstrap.Entity.Transformer;

const LiteMolPlugin = ({ plugin, setPlugin }) => {
  const pluginRef = useRef(null);
  useEffect(() => {
    const _plugin = LiteMol.Plugin.create({
      target: pluginRef.current,
      layoutState: {
        hideControls: true,
        collapsedControlsLayout:
          LiteMol.Bootstrap.Components.CollapsedControlsLayout.Landscape,
      },
      viewportBackground: '#F1F1F1',
    });
    setPlugin(_plugin);
  }, []);

  return (
    <div
      ref={pluginRef}
      style={{
        position: 'absolute',
        top: 0,
        right: 0,
        left: 0,
        bottom: 0,
      }}
    />
  );
};

const LiteMolTab = () => {
  const [plugin, setPlugin] = useState(null);

  useEffect(() => {
    if (!plugin) {
      return;
    }
    let action = plugin.createTransform();
    const pdbId = '6zgg';

    action
      .add(plugin.root, Transformer.Data.Download, {
        url: `https://www.ebi.ac.uk/pdbe/static/entry/${pdbId}_updated.cif`,
        type: 'String',
        id: pdbId,
      })
      .then(Transformer.Data.ParseCif, { id: pdbId }, { isBinding: true })
      .then(Transformer.Molecule.CreateFromMmCif, { blockIndex: 0 })
      .then(Transformer.Molecule.CreateModel, { modelIndex: 0 })
      .then(Transformer.Molecule.CreateMacromoleculeVisual, {
        polymer: true,
        het: true,
        water: false,
      });

    plugin.applyTransform(action);
  }, [plugin]);

  return <LiteMolPlugin plugin={plugin} setPlugin={setPlugin} />;
};

export default LiteMolTab;
