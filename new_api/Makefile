ifndef CC
	CC=gcc
endif

all: tools benchmarks dev

tools: htslib
	cd tools; make HTSLIB=`readlink -f $(HTSLIB)`
benchmarks: htslib
	cd benchmarks; make HTSLIB=`readlink -f $(HTSLIB)`
dev: htslib
	cd dev; make HTSLIB=`readlink -f $(HTSLIB)`

htslib:
	if [[ '$(HTSLIB)' == '' ]]; \
	then echo "Error: Please pass HTSLIB=... with path to htslib"; exit 1; fi

clean:
	rm -rf bin/*
	cd benchmarks; make clean
	cd dev; make clean

.PHONY: all clean tools benchmarks dev htslib
