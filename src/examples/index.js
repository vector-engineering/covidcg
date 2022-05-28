import { config } from '../config';
import rsv_examples from './examples.rsv';
import sars2_examples from './examples.sars2';

export function getExampleItems({ selectTree }) {
  switch (config.virus) {
    case 'rsv':
      return rsv_examples({ selectTree });
    case 'sars2':
      return sars2_examples({ selectTree });
  }
}
