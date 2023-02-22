# simple-demix\_check
Reimplementation of
[demix\_check](https://github.com/harry-thorpe/demix_check) using GNU
coreutils. This repository provides scripts to perform the steps to
set up the demix\_check reference and check reads binned with
[mGEMS](https://github.com/PROBIC/mGEMS) without requiring the python
or R dependencies in the original.

## Dependencies
Version numbers are for versions that have been tested. Older or newer
versions should also work.

- [seqtk](https://github.com/lh3/seqtk) v1.3
- [mash](https://github.com/marbl/Mash) v2.3
- [datamash](https://www.gnu.org/software/datamash/) v1.8

## Installation

If the dependencies are not installed on your system, run the supplied `install.sh` script to install them from
source (requires a C/C++ compiler). After running the script, you will
need to export the paths to these tools by
running the following command in the directory where simple-demix\_check is installed before running the scripts:
```
export PATH=$PATH:$(pwd)
```

## Usage
### Set up a set of reference sequences
```
nthreads=4
memmegas=1024M
tmpdir=tmp

setup_reference.sh --ref_info ref_info.tsv --threads $nthreads --tmpdir $tmpdir --bufzie $memmegas
```

### Check reads from mGEMS
``` abundances=cluster_abundances.txt nthreads=4 memmegas=1024M
tmpdir=tmp ref_dir=setup_reference_output forward=reads_1.fastq.gz
reverse=reads_2.fastq.gz

check_reads.sh --abundances $abundances --threads $nthreads --tmpdir $tmpdir --bufsize $memmegas --reference $ref_dir --fwd $forward --rev $reverse
```

## License
simple-demix_check is licensed under the [BSD-3-Clause
license](https://opensource.org/licenses/BSD-3-Clause). A copy of the
license is supplied with the project, or can alternatively be obtained
from
[https://opensource.org/licenses/BSD-3-Clause](https://opensource.org/licenses/BSD-3-Clause).
