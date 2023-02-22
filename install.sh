#!/bin/sh

set -euo pipefail

path=$(pwd)

## Compile datamash v1.8
wget http://ftp.gnu.org/gnu/datamash/datamash-1.8.tar.gz
tar -xzf datamash-1.8.tar.gz
cd datamash-1.8
./configure
make
datamash_path=$(pwd)"/"datamash

## Compile seqtk v1.3
cd ../
wget -O seqtk-v1.3.tar.gz https://github.com/lh3/seqtk/archive/refs/tags/v1.3.tar.gz
tar -xzf seqtk-v1.3.tar.gz
cd seqtk-1.3
make
seqtk_path=$(pwd)"/"seqtk

cd ../

## Downlaod mash v2.3 binary for Linux or MacOS depending on system
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    wget https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar
    tar -xf mash-Linux64-v2.3.tar
    mash_path=$(pwd)"/mash-Linux64-v2.3/mash"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    wget https://github.com/marbl/Mash/releases/download/v2.3/mash-OSX64-v2.3.tar
    tar -xf mash-OSX64-v2.3.tar
    mash_path=$(pwd)"/mash-OSX64-v2.3/mash"
fi

## Setup symlinks in root directory
cd ../
ln -s $mash_path $path"/"mash
ln -s $datamash_path $path"/"datamash
ln -s $seqtk_path $path"/"seqtk
