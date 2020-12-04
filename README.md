## COVID-19 CG (CoV Genetics)

**Preprint now up on bioRxiv: [https://www.biorxiv.org/content/10.1101/2020.09.23.310565v2](https://www.biorxiv.org/content/10.1101/2020.09.23.310565v2)**

Table of Contents

- [COVID-19 CG (CoV Genetics)](#covid-19-cg-cov-genetics)
- [Requirements](#requirements)
- [Installation](#installation)
  - [Python](#python)
    - [Data Requirements](#data-requirements)
    - [Data Package](#data-package)
  - [Javascript](#javascript)
    - [macOS](#macos)
    - [Linux](#linux)
- [Analysis Pipeline](#analysis-pipeline)
- [About the project](#about-the-project)
- [Data enabling COVID CG](#data-enabling-covid-cg)
- [Citing COVID CG](#citing-covid-cg)
  - [License](#license)
  - [Contributing](#contributing)


## Requirements

- Get a conda distribution of python, we recommend [miniconda3](https://docs.conda.io/en/latest/miniconda.html).
- curl
- node.js > 8.0.0
- npm

## Installation

`git clone https://github.com/vector-engineering/covidcg.git`

Install dependencies

```sh
conda env create -n covidcg -f environment.yml
```

### Javascript

This app was built from the [react-slingshot](https://github.com/coryhouse/react-slingshot) example app.

1. **Install [Node 8.0.0 or greater](https://nodejs.org)**

   Need to run multiple versions of Node? Use [nvm](https://github.com/creationix/nvm).

2. **Install [Git](https://git-scm.com/downloads)**.

3. **[Disable safe write in your editor](https://webpack.js.org/guides/development/#adjusting-your-text-editor)** to assure hot reloading works properly.

4. Complete the steps below for your operating system:

   #### macOS

   - Install [watchman](https://facebook.github.io/watchman/) via `brew install watchman` to avoid [this issue](https://github.com/facebook/create-react-app/issues/871) which occurs if your macOS has no appropriate file watching service installed.

   #### Linux

   - Run this to [increase the limit](http://stackoverflow.com/questions/16748737/grunt-watch-error-waiting-fatal-error-watch-enospc) on the number of files Linux will watch. [Here's why](https://github.com/coryhouse/react-slingshot/issues/6).

     `echo fs.inotify.max_user_watches=524288 | sudo tee -a /etc/sysctl.conf && sudo sysctl -p`.

5. **Install NPM packages**

   `npm install`

6. **Run the example app**

   `npm start -s`

   This will run the automated build process, start up a webserver, and open the application in your default browser. When doing development with this kit, this command will continue watching all your files. Every time you hit save the code is rebuilt, linting runs, and tests run automatically. Note: The -s flag is optional. It enables silent mode which suppresses unnecessary messages during the build.

---

## Analysis Pipeline

Data analysis is run with [Snakemake](https://snakemake.readthedocs.io/en/stable/), Python scripts, and bioinformatics tools such as `bowtie2`. Please ensure that the conda environment is configured correctly (See [Python](#Python)) and that all [data files](#Data-Requirements) are present and linked correctly to the `data/` folder.

Preview (dry-run) the pipeline by running:

```
snakemake -n
```

and run the pipeline with:

```
snakemake
```

**NOTE**: `bowtie2` usually uses anywhere from 8 – 10 GB of RAM per CPU during the alignment step. If the pipeline includes the alignment step, then only use as many cores as you have RAM / 10. i.e., if your machine has 128 GB RAM, then you can run at most 128 / 10 ~= 12 cores.

---

## About the project

This project is developed by the [Vector Engineering Lab](https://vector.engineering/):

- Albert Chen (Broad Institute)
- [Kevin Altschuler](https://www.linkedin.com/in/kevinaltschuler/)
- Shing Hei Zhan, PhD (University of British Columbia)
- Alina Yujia Chan, PhD (Broad Institute)
- Ben Deverman, PhD (Broad Institute)

The manuscript for this project is currently being prepared.

Contact the authors by email: [covidcg@broadinstitute.org](mailto:covidcg@broadinstitute.org)

Python/snakemake scripts were run and tested on MacOS 10.15.4 (8 threads, 16 GB RAM), Google Cloud Debian 10 (buster), (64 threads, 412 GB RAM), and Windows 10/Ubuntu 20.04 via. WSL2 (48 threads, 128 GB RAM)

## Data enabling COVID CG

We are extremely grateful to the [GISAID Initiative](https://www.gisaid.org/) and all its data contributors, i.e. the Authors from the Originating laboratories responsible for obtaining the speciments and the Submitting laboratories where genetic sequence data were generated and shared via the GISAID Initiative, on which this research is based.

Elbe, S., and Buckland-Merrett, G. (2017) Data, disease and diplomacy: GISAID’s innovative contribution to global health. _Global Challenges_, 1:33-46. DOI:[10.1002/gch2.1018](https://doi.org/10.1002/gch2.1018) PMCID: [31565258](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6607375/)

## Citing COVID CG

Users are encouraged to share, download, and further analyze data from this site. Plots can be downloaded as PNG or SVG files, and the data powering the plots and tables can be downloaded as well. Please attribute any data/images to [covidcg.org](https://covidcg.org/).

Note: When using results from these analyses in your manuscript, ensure that you acknowledge the contributors of data, i.e. _We gratefully acknowledge all the Authors from the Originating laboratories responsible for obtaining the speciments and the Submitting laboratories where genetic sequence data were generated and shared via the GISAID Initiative, on which this research is based_.

and cite the following reference(s):

Shu, Y., McCauley, J. (2017) GISAID: Global initiative on sharing all influenza data – from vision to reality. _EuroSurveillance_, 22(13) DOI:[10.2807/1560-7917.ES.2017.22.13.30494](https://doi.org/10.2807/1560-7917.ES.2017.22.13.30494) PMCID: [PMC5388101](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5388101/)

### License

COVID-19 CG is distributed by an [MIT license](https://github.com/vector-engineering/covidcg/blob/master/LICENSE).

### Contributing

Please feel free to contribute to this project by opening an issue or pull request in the [GitHub repository](https://github.com/vector-engineering/covidcg).
