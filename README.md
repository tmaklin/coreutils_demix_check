# simple-demix\_check
Reimplementation of
[demix\_check](https://github.com/harry-thorpe/demix_check) using GNU
coreutils. This repository provides scripts to perform the steps to
set up the demix\_check reference without using python or R.

## Dependencies
- seqtk v1.3
- mash v2.3
- GNU Parallel
- GNU datamash
- awk
- zcat
- pigz

## Usage
```
nthreads=4
setup_reference.sh ref_info.tsv $nthreads
```

## License
simple-demix_check is licensed under the [BSD-3-Clause
license](https://opensource.org/licenses/BSD-3-Clause). A copy of the
license is supplied with the project, or can alternatively be obtained
from
[https://opensource.org/licenses/BSD-3-Clause](https://opensource.org/licenses/BSD-3-Clause).
