# dbOTUc
[![Build Status](https://travis-ci.org/scwatts/otudistclust.svg?branch=master)](https://travis-ci.org/scwatts/otudistclust)
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

Efficient and parallelised implementation of [dbOTU3](https://github.com/swo/dbotu3).


## Table of contents
* [Table of contents](#table-of-contents)
* [Requirements](#requirements)
* [Installation](#installation)
* [Usage](#usage)
* [License](#license)


## Requirements
There are no requirements for using the pre-compiled static binaries on 64-bit linux distributions. Otherwise, there are several libraries which are required for building and running. For further information, see [Compiling from source](#compiling-from-source).


## Installing
`dbOTUc` can be installed via pre-compiled binaries or from source.


### GNU/Linux
For most 64-bit linux distributions (e.g. Ubuntu, Debian, RedHat, etc) the easiest way to obtain `dbOTUc` is via statically compiled binaries on the releases page. These binaries can be downloaded and run immediately without any setup as they have no dependencies.


### Compiling from source
Compiling from source requires these libraries and software:
```
C++11 (gcc-4.9.0+, clang-4.9.0+, etc)
OpenMP 4.0+
GNU Scientific Library 2.1+
GNU getopt
GNU make
GNU autoconf
```

After meeting the above requirements, cmompiling and installing `dbOTUc` from source can be done by:
```bash
git clone https://github.com/scwatts/otudistclust.git
cd otudistclust
./autogen.sh
./configure --prefix=/usr/
make
make install
```
Once completed, the `dbOTUc` executables can be run from the command line.


## Usage
Simple command invocation:
```bash
dbotuc -c otu_table.tsv -f otu_representative.fasta -o otu_table_clustered.tsv -m cluster_members.tsv -t 4
```

## License
[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)
