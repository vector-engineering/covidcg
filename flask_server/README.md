## Install Instructions

Installation/usage of the server is only tested using `conda` environments, but should work with a python virtual environment (`virtualenv`). An installation outside of an isolated environment is not recommended.

Install the dependencies:

```bash
pip install -r flask_server/requirements.txt
```

## Running

The server can only be run from the project root. `flask_server/serve.sh` sets all of the necessary environment variables to run the server in the development environment.

In addition, the `CONFIGFILE` environment variable needs to be set prior to running the server.

For example:

```bash
cd covidcg # project root
CONFIGFILE=config/config_genbank.yaml flask_server/serve.sh
```

