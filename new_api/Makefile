ifdef HTSLIB
	MARGS=HTSLIB=$(HTSLIB)
endif

all: tools benchmarks dev

tools:
	cd tools; make $(MARGS)
benchmarks:
	cd benchmarks; make $(MARGS)
dev:
	cd dev; make $(MARGS)

clean:
	rm -rf bin/*
	cd benchmarks; make clean
	cd dev; make clean

.PHONY: all clean tools benchmarks dev
