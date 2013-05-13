ifndef CC
	CC=gcc
endif

# Use bash as shell
SHELL := /bin/bash

CFLAGS=-Wall -Wextra -I .. -I $(HTSLIB)/htslib/ -L $(HTSLIB)/htslib/
LINKING=-lhts -lpthread -lz

ifdef DEBUG
	OPT = -O0 --debug -g
else
	OPT = -O3
endif

all: seqcmdline

seqcmdline: seq_cmdline.c ../seq_file.h ../stream_buffer.h | htslib
	$(CC) $(OPT) $(CFLAGS) -o seqcmdline seq_cmdline.c $(LINKING)

htslib:
	@if [[ '$(HTSLIB)' == '' ]]; \
	then echo "Error: Please pass HTSLIB=... with path to htslib"; exit 1; fi

clean:
	rm -rf seqcmdline *.dSYM *.greg

.PHONY: all clean htslib
