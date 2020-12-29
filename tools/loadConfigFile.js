import fs from 'fs';
import yaml from 'js-yaml';

let _configfile = {};
try {
  _configfile = yaml.safeLoad(fs.readFileSync(process.env.CONFIGFILE, 'utf8'));
  console.log('LOADED CONFIGURATION FILE');
  console.log(_configfile);
} catch (e) {
  console.error(e);
}

export const configfile = _configfile;
