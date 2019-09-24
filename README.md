# The shared neoantigen landscape of MSI cancers reflects immunoediting during tumor evolution

This repository contains all scripts used in the preprint "The shared neoantigen landscape of MSI cancers reflects immunoediting during tumor evolution". This includes data anonymization, processing and generation of figures.

To re-generate data and plots, run `run.sh` from this directory.

## Analysis setup

To run the script for rerunning our analysis on processed data, first perform the following steps.
In addition to standard GNU/Linux tools the following are necessary for use.

### Install dependencies

Set up a conda environment:

conda env create --name neomsi --file=env.yml
conda activate neomsi

Install NetMHCpan (instructions: http://www.cbs.dtu.dk/services/doc/netMHCpan-4.0.readme), and add `netMHCpan` to your path:

```bash
export PATH=/path/to/netmhcpan/bin/:$PATH
export NETMHCpan=/path/to/netmhcpan/<Linux/Darwin>/
```

once installation is done, you can run `run.sh` to rerun analysis and figure generation.

## ReFrame setup

Change directory to `ReFrame`. To run ReFrame on example data, run the following:

```bash
r -f reframe.R
```

This will process all `.xlsx` files containing raw fragment length data in the `in` directory and deposit processed
allele frequency data to the `out` folder. Example input and output files are provided in the respective folders.
