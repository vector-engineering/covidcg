![](https://covidcg.org/cg_logo_v13.png)

## COVID-19 CG (CoV Genetics)

**Article now up at eLife: [https://doi.org/10.7554/eLife.63409](https://doi.org/10.7554/eLife.63409)**

Table of Contents

- [COVID-19 CG (CoV Genetics)](#covid-19-cg-cov-genetics)
- [Data enabling COVID CG](#data-enabling-covid-cg)
- [Installation](#installation)
  - [Dependency changes](#dependency-changes)
  - [Database refresh](#database-refresh)
- [Per-service installation](#per-service-installation)
  - [Javascript](#javascript)
  - [PostgreSQL](#postgresql)
  - [Flask Server](#flask-server)
- [Analysis Pipeline](#analysis-pipeline)
  - [Pipeline Installation](#pipeline-installation)
  - [Ingestion](#ingestion)
  - [Main Analysis](#main-analysis)
- [About the project](#about-the-project)
- [Citing COVID CG](#citing-covid-cg)
  - [License](#license)
  - [Contributing](#contributing)

## Data enabling COVID CG

We are extremely grateful to the [GISAID Initiative](https://www.gisaid.org/) and all its data contributors, i.e. the Authors from the Originating laboratories responsible for obtaining the speciments and the Submitting laboratories where genetic sequence data were generated and shared via the GISAID Initiative, on which this research is based.

Elbe, S., and Buckland-Merrett, G. (2017) Data, disease and diplomacy: GISAID’s innovative contribution to global health. _Global Challenges_, 1:33-46. DOI:[10.1002/gch2.1018](https://doi.org/10.1002/gch2.1018) PMCID: [31565258](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6607375/)

Note: When using results from these analyses in your manuscript, ensure that you acknowledge the contributors of data, i.e. _We gratefully acknowledge all the Authors from the Originating laboratories responsible for obtaining the speciments and the Submitting laboratories where genetic sequence data were generated and shared via the GISAID Initiative, on which this research is based_.

## Installation

The COVID-19 CG website comprises of 3 services (PostgreSQL database, Flask server, React frontend). These can be run separately (see detailed instructions at [per-service installation](#per-service-installation)) but we recommend using Docker to manage these services.

The analysis pipeline for processing raw SARS-CoV-2 genomes is a separate install, and described below in [Analysis Pipeline](#analysis-pipeline)

- Install [Docker](https://docs.docker.com/get-docker/)
- Clone this repository: `git clone https://github.com/vector-engineering/covidcg.git`

```bash
$ cd covidcg
$ docker-compose build # Build containers
                       # (Re-builds only necessary if packages or
                       # dependencies have changed)
$ docker-compose up -d # Run all services
$ docker-compose down # Shut down all services when finished
```

The default deployment (`docker-compose.yml`) will run all 3 sites at the same time (sars2, rsv, and flu). For virus-specific sites, see `docker-compose.sars2.yml`, etc. Run a specific deployment with:

```bash
docker compose -f docker-compose.sars2.yml build
docker compose -f docker-compose.sars2.yml up -d
...
```

**NOTE**: When starting from a fresh database, the server will automatically seed the database with data from the `example_data_genbank` folder. Data provided with the repository is in raw/gzipped form and needs to be unarchived and processed before the data can be loaded into the database. Please see the [Analysis Pipeline](#analysis-pipeline) section for instructions on processing this data.

### Dependency changes

If the dependencies for the JS change (i.e., a change in `package.json`), then you can rebuild the `cg-frontend` container with:

```bash
$ docker-compose down
$ docker-compose build --no-cache cg-frontend
$ docker-compose up
```

A rebuild will also need to be run if the toolchains change (`webpack*.js` or anything in `tools/`)

For files outside of `src`, i.e., in `config/` or in `static_data/`, the container will need to be restarted but not rebuilt.

For dependency changes for the server (i.e., changes in `requirements.txt`)

```bash
$ docker-compose down
$ docker-compose build --no-cache cg-server
$ docker-compose up
```

### Database refresh

To erase the local development database, delete the postgres docker volume with:

```bash
$ docker-compose down -v # -v will delete the volume
$ docker-compose up
```

## Per-service installation

We recommend developing with Docker and `docker-compose`. More details on the installation for each service can be found in their respective `Dockerfile`s in the `services/` folder, and in the `docker-compose.yml` file. Running each service separately is not recommended and not tested on our end. Since we are not actively testing per-service installations, please submit a GitHub issue if you run into any problems during installation or running.

First, clone this repository: `git clone https://github.com/vector-engineering/covidcg.git`

### Javascript

Requirements:

- curl
- node.js > 8.0.0
- npm

This app was built from the [react-slingshot](https://github.com/coryhouse/react-slingshot) example app.

1. **Install [Node 8.0.0 or greater](https://nodejs.org)**

   Need to run multiple versions of Node? Use [nvm](https://github.com/creationix/nvm).

2. **Install [Git](https://git-scm.com/downloads)**.

3. **[Disable safe write in your editor](https://webpack.js.org/guides/development/#adjusting-your-text-editor)** to assure hot reloading works properly.

4. Complete the steps below for your operating system:

   macOS

   - Install [watchman](https://facebook.github.io/watchman/) via `brew install watchman` to avoid [this issue](https://github.com/facebook/create-react-app/issues/871) which occurs if your macOS has no appropriate file watching service installed.

   Linux

   - Run this to [increase the limit](http://stackoverflow.com/questions/16748737/grunt-watch-error-waiting-fatal-error-watch-enospc) on the number of files Linux will watch. [Here's why](https://github.com/coryhouse/react-slingshot/issues/6).

     `echo fs.inotify.max_user_watches=524288 | sudo tee -a /etc/sysctl.conf && sudo sysctl -p`.

5. **Install NPM packages**

   `npm install`

6. **Run the app**

   `CONFIGFILE=config/config_genbank.yaml npm start -s`

   This will run the automated build process, start up a webserver, and open the application in your default browser. When doing development with this kit, this command will continue watching all your files. Every time you hit save the code is rebuilt, linting runs, and tests run automatically. Note: The -s flag is optional. It enables silent mode which suppresses unnecessary messages during the build.

### PostgreSQL

This development environment was tested with PostgreSQL 12

Please provide DB connection information to the Flask server with the following environment variables:

- POSTGRES_USER
- POSTGRES_PASSWORD
- POSTGRES_DB
- POSTGRES_HOST
- POSTGRES_PORT
- POSTGRES_MAX_CONN (the maximum number of connections for the Postgres connection pool)

### Flask Server

Requirements:

- Python3 (Python >= 3.8) with virtual environments. We recommend conda via. [miniconda3](https://docs.conda.io/en/latest/miniconda.html), but python3 with `virtualenv` or any other virtual environment provider should also work fine

Install dependencies:

```bash
$ cd services/server
$ pip install -r requirements.txt
```

Run server:

```bash
$ cd services/server
$ CONFIGFILE=../../config/config_genbank.yaml ./serve.sh # Run Flask server in development mode, with GenBank settings
                                                   # Optionally, edit the serve.sh script to set the config file
```

---

## Analysis Pipeline

Data analysis is run with [Snakemake](https://snakemake.readthedocs.io/en/stable/), Python scripts, and bioinformatics tools such as `bowtie2`. Please ensure that the conda environment is configured correctly (See [Pipeline Installation](#Pipeline-Installation)).

Data analysis is broken up into two snakemake pipelines: 1) ingestion and 2) main. The ingestion pipeline downloads, chunks, and prepares metadata for the main analysis, and the main pipeline analyzes sequences, extracts mutations, and compiles data for display in the web application.

Configuration of the pipeline is defined in the `config/config_[workflow].yaml` files.

### Pipeline Installation

1. Clone this repository: `git clone https://github.com/vector-engineering/covidcg.git`
2. Install [miniconda3](https://docs.conda.io/en/latest/miniconda.html)
3. Create conda environment:

```bash
$ conda config --add channels bioconda # Add package download locations
$ conda config --add channels conda-forge
$ conda env create -f environment.yml
```

For OSX M1 chips, use the alternative environment `environment_osx-arm64.yaml`. Some additional source compilation steps are required as not all ARM64 binaries are available on conda.

### Ingestion

Currently available ingest workflows are:

SARS2:

- `workflow_sars2_gisaid_ingest`
- `workflow_sars2_genbank_ingest`
- `workflow_sars2_custom_ingest`

RSV:

- `workflow_rsv_genbank_ingest`
- `workflow_rsv_custom_ingest`

Flu:

- `workflow_flu_genbank_ingest`
- `workflow_flu_custom_ingest`

**NOTE: While GISAID ingestion pipelines are provided as open-source, it is intended only for internal use**.

GenBank ingest pipelines are designed to automatically download and process data from their respective data source.

"Custom" ingest pipelines can be used for analyzing and visualizing in-house data. More details are available in README files within each ingestion pipeline's folder. Each ingestion workflow is parametrized by its own config file. i.e., `config/config_sars2_genbank.yaml` for the SARS-CoV-2 GenBank workflow.

For example, you can run the SARS-CoV-2 GenBank ingestion pipeline with:

```bash
$ cd workflow_sars2_genbank_ingest
$ snakemake --use-conda # Conda required specifically for SARS2 GenBank ingest in order to run Pangolin lineage assignments
```

### Main Analysis

The main data analysis pipeline is located in `workflow_main`. It requires data, in a data folder, from the ingestion pipeline. The data folder is defined in the `config/config_[workflow].yaml` file. The path to the config file is required for the main workflow, as it needs to know what kind of data to expect (as described in the config files).

For example, if you ingested data from GenBank, run the main analysis pipeline with:

```bash
cd workflow_main
snakemake --configfile ../config/config_sars2_genbank_dev.yaml
```

This pipeline will align sequences to the reference sequence with `minimap2`, extract mutations on both the NT and AA level, and combine all metadata and mutation information data. The output data can be uploaded to a PostgreSQL database with `workflow_main/scripts/push_to_database.py`. Or, you can use the output files directly for your own analyses.

### Example data

Example data from GenBank is provided for all viruses, and is located in gzipped tarballs inside the `example_data_genbank` folder. Data for some viruses is truncated by submission date in order to lighten data load and speed up development on smaller machines.

To extract the data:

```bash
$ cd example_data_genbank
$ tar -xzf sars2.tar.gz
$ tar -xzf rsv.tar.gz
$ tar -xzf flu.tar.gz
```

These tarballs contain only raw sequences and metadata, and mimic the output from their respective ingest pipelines. Once the files are extracted, run the main analysis workflow described above.

---

## About the project

This project is developed by the [Vector Engineering Lab](https://vector.engineering/):

- Albert Tian Chen (Broad Institute)
- [Kevin Altschuler](https://www.linkedin.com/in/kevinaltschuler/)
- Shing Hei Zhan, PhD (University of British Columbia)
- Alina Yujia Chan, PhD (Broad Institute)
- Ben Deverman, PhD (Broad Institute)

Contact the authors by email: [covidcg@broadinstitute.org](mailto:covidcg@broadinstitute.org)

Python/snakemake scripts were run and tested on MacOS 10.15.4 (8 threads, 16 GB RAM), Google Cloud Debian 10 (buster), (64 threads, 412 GB RAM), and Windows 10/Ubuntu 20.04 via. WSL2 (48 threads, 128 GB RAM)

## Citing COVID CG

Users are encouraged to share, download, and further analyze data from this site. Plots can be downloaded as PNG or SVG files, and the data powering the plots and tables can be downloaded as well. Please attribute any data/images to [covidcg.org](https://covidcg.org/), or cite our manuscript:

Chen AT, Altschuler K, Zhan SH, Chan YA, Deverman BE. COVID-19 CG enables SARS-CoV-2 mutation and lineage tracking by locations and dates of interest. _eLife_ (2021), doi: [https://doi.org/10.7554/eLife.63409](https://doi.org/10.7554/eLife.63409)

Note: When using results from these analyses in your manuscript, ensure that you acknowledge the contributors of data, i.e. _We gratefully acknowledge all the Authors from the Originating laboratories responsible for obtaining the speciments and the Submitting laboratories where genetic sequence data were generated and shared via the GISAID Initiative, on which this research is based_.

and cite the following reference(s):

Shu, Y., McCauley, J. (2017) GISAID: Global initiative on sharing all influenza data – from vision to reality. _EuroSurveillance_, 22(13) DOI:[10.2807/1560-7917.ES.2017.22.13.30494](https://doi.org/10.2807/1560-7917.ES.2017.22.13.30494) PMCID: [PMC5388101](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5388101/)

### License

COVID-19 CG is distributed by an [MIT license](https://github.com/vector-engineering/covidcg/blob/master/LICENSE).

### Contributing

Please feel free to contribute to this project by opening an issue or pull request in the [GitHub repository](https://github.com/vector-engineering/covidcg).
