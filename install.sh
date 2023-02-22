wget http://ftp.gnu.org/gnu/datamash/datamash-1.8.tar.gz
tar -xzf datamash-1.8.tar.gz
cd datamash-1.8
./configure
make

cd ../
wget -O seqtk-v1.3.tar.gz https://github.com/lh3/seqtk/archive/refs/tags/v1.3.tar.gz
tar -xzf seqtk-v1.3.tar.gz
cd seqtk-1.3
make

cd ../
wget https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar
tar -xf mash-Linux64-v2.3.tar

