#!/bin/bash

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

## Compile sketchlib from source
git clone https://github.com/bacpop/sketchlib.rust
cd sketchlib.rust
git checkout 37a5879
cargo build --release
sketchlib_path=$(pwd)"/target/release/sketchlib"

## Setup symlinks in root directory
cd ../
ln -s $sketchlib_path $path"/"sketchlib
ln -s $datamash_path $path"/"datamash
ln -s $seqtk_path $path"/"seqtk
