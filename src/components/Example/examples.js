import { config } from '../../config';
import sars2 from './examples.sars2';

export function getExampleItems({ selectTree }) {
  switch (config.virus) {
    case 'rsv':
      return [];
    case 'sars2':
      return sars2(selectTree);
    default:
      return import('./examples.sars2');
  }
}
