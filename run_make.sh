#! /bin/sh
mkdir -p lib
export EIGEN3_INCLUDE_DIR=$(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix'])")/include/eigen3
make
cp lib/* ../sharpy/lib/
