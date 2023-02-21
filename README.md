# simple-demix\_check
Reimplementation of
[demix\_check](https://github.com/harry-thorpe/demix_check) using GNU
coreutils. This repository provides scripts to perform the steps to
set up the demix\_check reference and check reads binned with
[mGEMS](https://github.com/PROBIC/mGEMS) without requiring the python
or R dependencies in the original.

## Dependencies
- seqtk v1.3
- mash v2.3
- GNU Parallel
- GNU datamash
- awk
- zcat
- pigz

## Usage
### Set up a set of reference sequences
```
nthreads=4
memmegas=1024
tmpdir=tmp

setup_reference.sh ref_info.tsv $nthreads $tmpdir $memmegas"M"
```

### Check reads from mGEMS
```
abundances=cluster_abundances.txt
nthreads=4
memmegas=1024
tmpdir=tmp
ref_dir=setup_reference_output
forward=reads_1.fastq.gz
reverse=reads_2.fastq.gz

check_reads.sh $abundances $nthreads $tmpdir $memmegas"M" $ref_dir $forward $reverse
```

## License
simple-demix_check is licensed under the [BSD-3-Clause
license](https://opensource.org/licenses/BSD-3-Clause). A copy of the
license is supplied with the project, or can alternatively be obtained
from
[https://opensource.org/licenses/BSD-3-Clause](https://opensource.org/licenses/BSD-3-Clause).
