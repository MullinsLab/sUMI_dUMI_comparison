
This is an early version of PORPIDpipeline that uses the sample ID and UMI of the sUMI
primer for initial demultiplexing, QC and template identification. The "likely-real" families
from UMI1 are then passed on for UMI2 identification and grouping using the dUMI primer sequence.
The pipeline provides information on UMI combinations and frequency in the tagged directory as well as
consensus sequences derived from all reads in a UMI1 collection (sUMI) or only the most 
prevalent UMI1-UMI2 combination (dUMI).

This pipeline was created to provide matched sUMI and dUMI datasets for the analysis in
the manuscript "Optimized SMRT-UMI protocol produces highly accurate sequence datasets 
from diverse populations – application to HIV-1 quasispecies" Westfall et al.

PORPIDpipeline
https://github.com/MurrellGroup/PORPIDpipeline.git


## Quick start

### Dependencies (on an ubuntu machine)

- first update all apps
   - `apt update`
   - `apt upgrade`
- Snakemake
   - `apt-get install -y snakemake`
- mafft
   - `apt-get install -y mafft`
- fasttree
   - `apt-get install -y fasttree`
- python3 packages
  - `apt-get install python3-pandas`
  - `apt-get install python3-seaborn`


### Julia version 1.7

Download and unpack the latest Julia (we recommend version 1.7.1) from: 

[https://julialang.org/downloads/](https://julialang.org/downloads/)

Make sure you can enter the julia REPL from the command line, on an ubuntu machine you would do:

```bash
# move the julia system to a lib directory
mv julia-1.7.1 /usr/lib/julia-1.7.1
# make julia v1.7.1 executable from the command line
ln -s /usr/lib/julia-1.7.1/bin/julia /usr/local/bin/julia
# check that you can enter the julia REPL
julia --version
```

### cloning the PorpidPostproc repository

Now that the dependencies are setup we clone the PorpidPostproc repository

```bash
cd ~
git clone -b https://github.com/MullinsLab/sUMI_dUMI_comparison.git
```

### setting up the Julia package environment

then navigate to the `sUMI_dUMI_comparison` project folder and start the Julia REPL. 
Enter the package manager using `]` and then enter

```julia
activate .
instantiate
precompile
```

This will activate, install, and precompile the `julia` environment specified by the 
`Project.toml` and `Manifest.toml` files. The `precompile` command
above is not strictly needed but is useful if there are issues with installing
the `julia` packages listed in `Project.toml`

Next, add the following text to your Julia startup file (typically at `~/.julia/config/startup.jl`; 
you may need to create the directory if not present, `mkdir -p ~/.julia/config`).

```julia
using Pkg
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
```

This will activate the local environment at Julia startup.


### Configuration

For this project, libraries are expected to contain a forward (sUMI) and reverse (dUMI) UMI primer. 
The config file is expected to follow this formatting scheme (see example_config.yaml):

```
dataset001:
  template001:
    Forward_Primer_2ndRd_Sequence: GATCTACACTCTTTCCCTACACGAC
    Reverse_Primer_2ndRd_Sequence: AGATGTGTATAAGAGACAGCCGCTCCGTCCGACGACTCACTATA
    sUMI_Primer_Sequence: actataNNNNNNNNTTAGAGACATCCCCAGAGCTGTTAG
    dUMI_Primer_Sequence: NNNNNNNNATAGCTGTCTTTTATC
```


The PCR and UMI primer sequences provided will be used for demultiplexing and will be trimmed
from the final sequences. **Forward_Primer_2ndRd_Sequence** and **Reverse_Primer_2ndRd_Sequence**
are the 2nd rd PCR primers. 

CCS .fastq files should be placed in the `fastq/` subdirectory and named 
according to the the dataset name used in the `example_config.yaml` file, ie, `Dataset1.fastq`
for the example dataset.


### Preview and Execution

Preview jobs with Snakemake and run with {n} cores.

```bash
#preview jobs
snakemake -np

#run
snakemake -j{n}
```

For more info on Snakemake, see https://snakemake.readthedocs.io/en/stable/

## Conda setup

Some (without root access) may prefer to setup PorpidPostproc in a **conda** environment.

To accomplish this, first install `anaconda` locally. (the install script allows you to choose
the location for anaconda, by default `/home/user` but choose something else if
you want something accessable to a group of users)

```bash
curl –O https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh > Anaconda3-2021.05-Linux-x86_64.sh
bash Anaconda3-2021.05-Linux-x86_64.sh
```

then log out and log in again and check that you are in the `base` environment.

`conda` is very slow, so we suggest installing `mamba` in the conda `base` environment:

```bash
conda install -n base -c conda-forge mamba
```
clone the PorpidPostproc repository

```bash
cd ~  # or some other directory used for your anaconda installation
git clone -b h705mod1 https://gitlab.com/hugh.murrell/porpidpostproc.git
```

and then all the PorpidPostproc dependencies including `julia` version `1.7.1`
( as listed in the PorpidPostproc conda environment spec in `environment.yaml`),
can be installed in a `conda` environment via `mamba` using the commands:

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
mamba env create --file environment.yaml
```

Note that if you did use *some other directory* than your home directory for
installing the PorpidPostproc repository then you have to inform Julia where
your packages are stored by placing the following command in your `.bashrc`
file:

```bash
# set path to .julia files
export JULIA_DEPOT_PATH="/some/other/directory/.julia"
```

to complete the setup, activate the new PorpidPostproc conda environment, 

```bash
conda activate PorpidPostproc
```

and continue with the `julia` package environment setup as outlined above in the *quick start* section.