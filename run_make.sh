#! /bin/sh
source activate sharpy
mkdir -p lib
make
cp lib/* ../sharpy/lib/
