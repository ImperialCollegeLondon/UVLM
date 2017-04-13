UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
EXT = so
endif
ifeq ($(UNAME), Darwin)
EXT = dylib
endif
default:
	mkdir -p lib
	$(MAKE) -C src/

clean:
	cd src; make clean
