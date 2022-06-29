import { config } from '../config';
import rsv_examples from './examples.rsv';
import sars2_examples from './examples.sars2';
import flu_examples from './examples.flu';

export function getExampleItems({ selectTree }) {
  switch (config.virus) {
    case 'rsv':
      return rsv_examples({ selectTree });
    case 'sars2':
      return sars2_examples({ selectTree });
    case 'flu':
      return flu_examples({ selectTree });
  }
}
