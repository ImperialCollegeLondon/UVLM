## LINUX G++ SUPPORT
FLAGS = -fPIC -O3 -march=native -std=c++14 -I$(EIGEN3_INCLUDE_DIR) -fomit-frame-pointer -ffast-math -fopenmp
LINKER_FLAGS = -shared -fPIC -fopenmp

## MACOS G++ SUPPORT
# FLAGS = -fPIC -g -std=c++14 -I$(EIGEN3_INCLUDE_DIR) -fomit-frame-pointer
# LINKER_FLAGS = -shared -fPIC

## INTEL C++ MKL SUPPORT
FLAGS = -fPIC -O3 -march=native -std=c++14 -I$(EIGEN3_INCLUDE_DIR) -fomit-frame-pointer -ffast-math -mkl=parallel
LINKER_FLAGS = -shared -fPIC -liomp5 -lpthread -ldl -lm 


UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
# CPP = icc
CPP = g++-5
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
