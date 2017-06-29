# FLAGS = -g3 -fPIC -O0 -Wall -std=c++14 -I/usr/local/eigen/
# FLAGS = -fPIC -O3 -march=native -std=c++14 -I/usr/local/eigen/# -fopenmp
FLAGS = -fPIC -O3 -march=native -std=c++14 -I/usr/local/eigen/ -fomit-frame-pointer# -fopenmp
LINKER_FLAGS = -shared -fPIC #-fopenmp
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
# FC = gfortran
CPP = g++-5
endif
ifeq ($(UNAME), Darwin)
# FC = ifort
CPP = g++-5
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
