
This is an early version of the porpid pipeline that uses the sample ID and UMI of the sUMI
primer for initial demultiplexing, QC and template identification. The "likely-real" families
from UMI1 are then passed on for UMI2 identification and grouping using the dUMI primer sequence.
The pipeline provides information on UMI combinations and frequency in the tagged directory as well as
consensus sequences derived from all reads in a UMI1 collection (sUMI) or only the most 
prevalent UMI1-UMI2 combination (dUMI).


### Setup

For this project, libraries are expected to contain a forward (sUMI) and reverse (dUMI) UMI primer. The config file is expected to follow this formatting scheme:

```
dataset001:
  template001:
    Forward_Primer_2ndRd_Sequence: GATCTACACTCTTTCCCTACACGAC
    Reverse_Primer_2ndRd_Sequence: AGATGTGTATAAGAGACAGCCGCTCCGTCCGACGACTCACTATA
    sUMI_Primer_Sequence: actataNNNNNNNNTTAGAGACATCCCCAGAGCTGTTAG
    dUMI_Primer_Sequence: NNNNNNNNATAGCTGTCTTTTATC
```

The file `config.yaml` contains an example. The supplied primer sequences are trimmed before UMI calling.

### Usage

Ensure all dependencies are installed, navigate to the project directory and run

```{bash}
snakemake --cores n
```

where "n" is the number of concurrent jobs you wish to run. For more on snakemake, look [here](https://snakemake.readthedocs.io/en/stable/).

### Dependencies

Julia >v1.0 must be installed on your system, which you can download from [here](https://julialang.org/downloads/). If they are not already present on your system, fire up the Julia interpreter, hit "]"  to enter the package manager and run the following cell to install the required Julia packages from the supplied Project.toml/Manifest.toml files:

```{julia}
activate .
instantiate
precompile
```

This pipeline also provides a conda environment to set up the other dependencies, including snakemake. Install conda and run

```
conda env create --file environment.yaml
conda activate porpid
```
