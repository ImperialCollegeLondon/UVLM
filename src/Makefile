## LINUX G++ SUPPORT
#FLAGS = -fPIC -O3 -march=native -std=c++14 -I$(EIGEN3_INCLUDE_DIR) -fomit-frame-pointer -ffast-math -fopenmp -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE
#LINKER_FLAGS = -L$(PREFIX)/lib -shared -fPIC -fopenmp -lblas -llapacke -Wl,-z,defs

## G++ DEBUG (VERY SLOW)
#FLAGS = -fPIC -O0 -g3 -march=native -std=c++14 -I$(EIGEN3_INCLUDE_DIR) 
#LINKER_FLAGS = -L$(EIGEN3_INCLUDE_DIR)/../../lib -shared -fPIC -lblas -llapacke -Wl,-z,defs

## LINUX G++ + MKL
FLAGS = -fPIC -O3 -march=native -std=c++14 -DEIGEN_USE_MKL_ALL -I$(PREFIX)/include -I$(EIGEN3_INCLUDE_DIR) -fomit-frame-pointer -ffast-math -fopenmp
LINKER_FLAGS = -shared -fPIC -fopenmp -L$(PREFIX)/lib -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl


UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
# CPP = icc
# CPP = g++-5
CPP = g++
endif
ifeq ($(UNAME), Darwin)
CPP = g++
# CPP = clang
endif


LINKER = ar cr
object_dir = obj/
lib_dir = ../lib/


ifeq ($(UNAME), Linux)
target = libuvlm.so
endif
ifeq ($(UNAME), Darwin)
target = libuvlm.dylib
endif

files = cpp_interface

objects = $(files:=.o)
sources = $(files:=.cpp)

include_dir = ../include

#LINKER_FLAGS += -l$(base_lib) -L$(libs_dir)

default: $(target)

$(target): $(objects)
	$(CPP) -o $(target) $(addprefix $(object_dir), $(objects)) -I$(include_dir) $(LINKER_FLAGS)
	mv $(target) $(lib_dir)$(target)

%.o: %.cpp
	$(CPP) $(FLAGS) -I$(include_dir) -c $<
	mkdir -p $(object_dir)
	mv $@ $(object_dir)$@

clean:
	rm -rf $(object_dir)

.PHONY: default clean
