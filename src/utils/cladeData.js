import cladeData from '../../processed_data/clade_data.json';
import _ from 'underscore';

export function loadCladeData() {
  return cladeData;
}

export function getReferenceSequence() {
  return _.findWhere(cladeData, { clade: 'root' })['seq'];
}

export function getCladesFromGene(gene) {
  let start_pos = gene.start;
  let end_pos = gene.end;

  let cladeData = loadCladeData();
  let validClades = [];
  let pos = -1;

  // Get clades whose positions fall within start -- end
  cladeData.forEach(clade => {
    // If any one of it's positions is in the range, then add it
    for(let i = 0; i < clade.pos.length; i++) {
      pos = clade.pos[i];
      if(pos >= start_pos && pos <= end_pos) {
        validClades.push(clade);
        return; // Onto the next clade
      }
    }
  });
  
  return validClades;
}