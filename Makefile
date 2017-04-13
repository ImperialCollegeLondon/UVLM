UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
EXT = so
endif
ifeq ($(UNAME), Darwin)
EXT = dylib
endif
default:
	$(MAKE) -C src/
	mkdir -p lib
	mv src/lib/*.$(EXT) lib/

clean:
	cd src; make clean
