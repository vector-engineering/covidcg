# COVID-UI
(working title)

## Installation

`git clone --recursive https://github.com/vector-engineering/covid-ui.git`

### Python

1. Get a conda distribution of python, we recommend [miniconda3](https://docs.conda.io/en/latest/miniconda.html).
   
2. Install dependencies
    ```sh
    conda env create -f environment.yml
    ```

3. Install pangolin
    ```sh
    cd pangolin
    pip install --editable .
    ```

### Javascript

This app was built from the [react-slingshot](https://github.com/coryhouse/react-slingshot) example app.

1. **Install [Node 8.0.0 or greater](https://nodejs.org)**

    Need to run multiple versions of Node? Use [nvm](https://github.com/creationix/nvm).

2. **Install [Git](https://git-scm.com/downloads)**.

3. **[Disable safe write in your editor](https://webpack.js.org/guides/development/#adjusting-your-text-editor)** to assure hot reloading works properly.

4. Complete the steps below for your operating system:

    ### macOS

    * Install [watchman](https://facebook.github.io/watchman/) via `brew install watchman` to avoid [this issue](https://github.com/facebook/create-react-app/issues/871) which occurs if your macOS has no appropriate file watching service installed.

    ### Linux

    * Run this to [increase the limit](http://stackoverflow.com/questions/16748737/grunt-watch-error-waiting-fatal-error-watch-enospc) on the number of files Linux will watch. [Here's why](https://github.com/coryhouse/react-slingshot/issues/6).

        `echo fs.inotify.max_user_watches=524288 | sudo tee -a /etc/sysctl.conf && sudo sysctl -p`.

5. **Run the example app**

    `npm start -s`

    This will run the automated build process, start up a webserver, and open the application in your default browser. When doing development with this kit, this command will continue watching all your files. Every time you hit save the code is rebuilt, linting runs, and tests run automatically. Note: The -s flag is optional. It enables silent mode which suppresses unnecessary messages during the build.

---

## About the project

This project was developed by ...

The paper for this project ...

Contact the authors by email: ...

Python scripts were run on MacOS 10.15.4 (8 threads, 16 GB RAM) and Google Cloud Debian 10 (buster), (64 threads, 240 GB RAM)

## Acknowledgements

This project is powered by many open-source software projects.

Find all acknowledgements at ...

### License

... is distributed by an [MIT license](https://github.com/vector-engineering/covid-ui/blob/master/LICENSE.txt).

### Contributing

Please feel free to contribute to this project by opening an issue or pull request in the [GitHub repository](https://github.com/vector-engineering/covid-ui).
