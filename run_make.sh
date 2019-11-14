#! /bin/sh
mkdir -p lib
export PREFIX=$(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix'])")
if [ "$PREFIX" = "None" ]; then
    echo "*** Please check that the python environment is active."
    echo "*** Run ``conda activate sharpy_env``."
    exit 1
fi
export EIGEN3_INCLUDE_DIR=$PREFIX/include/eigen3
export MKL_ROOT=$PREFIX
make
cp lib/* ../sharpy/lib/
