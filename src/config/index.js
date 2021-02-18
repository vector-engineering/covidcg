// The configuration YAML file (defined by the environment variable CONFIGFILE
// when the app build is launched) is injected in place of the string "CG_CONFIG"
// as a JSON object, via webpack.DefinePlugin
const _config = CG_CONFIG; // eslint-disable-line no-undef

export const config = Object.assign({}, _config);

const _hostname =
  process.env.NODE_ENV == 'development'
    ? config['dev_hostname']
    : config['prod_hostname']; // eslint-disable-line no-undef

export const hostname = _hostname;
export const init_endpoint = _hostname + config['init_endpoint'];
export const data_endpoint = _hostname + config['data_endpoint'];
