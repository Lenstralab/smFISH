# smFISH
Lenstra labs code for single molecule FISH as used by Gowthaman et el, Brouwer et al and Patel et al.
This script uses Python 3.

## Cloning the repository
Install git: https://git-scm.com/

    git clone https://github.com/Lenstralab/smFISH.git
    cd smFISH

# Version of code used for Patel et al.
    git checkout 011e50731237ca2fc8427099996a788244fe5ec9

# Version of code used for Brouwer et al.
    git checkout 28e6f801bbf95c88a18d09074adbd96ca12b430a

# Version of code used for Gowthaman et al.
    git checkout 36e73147e1c61a57e3ebfac595017d01a0553130

## Installation
If not done already:
- Install python (at least 3.7): https://www.python.org
- Install pip and git
- First install SimpleElastix: https://simpleelastix.readthedocs.io/GettingStarted.html

Then install the smFISH script (up to 5 minutes):

    pip install -e . --user

This will install the smfish package in your personal path in editable mode.

## Running the pipeline
From a terminal:

    FISH_pipeline /path/to/parameter_file.yml

Note that because of the pip step in the installation the FISH_pipeline script is in your path
and accessible from anywhere. So you don't have to specify the path. If this is not what you
want, you can still call the script by path, just notice it's in the subfolder 'smfish':

    /path/to/repo/smfish/FISH_pipeline.py /path/to/parameter_file.yml

Or from an ipython window:

    %run /path/to/repo/smfish/FISH_pipeline.py /path/to/parameter_file.yml

## Running the demo

    cd smFISH_analysis
    FISH_pipeline test_files/FISH_pipeline_test.yml

The script will go through the various stages of the analysis and normally finish in about 10 minutes.
It will make a folder 'demo_output' (defined in the parameter .yml file) which contains results in .tif, .txt and .pdf
formats.

## Testing
This script was tested with python 3.8 on Ubuntu 20.04
