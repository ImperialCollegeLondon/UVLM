#! /bin/sh

mkdir -p lib
make
cp lib/* ../sharpy/lib/
