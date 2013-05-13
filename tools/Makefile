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

all: ../bin/facat ../bin/fqcat ../bin/seqcat ../bin/seqstat \
	   ../bin/revcmp ../bin/seqrnd

../bin/facat: seq_cat.c ../seq_file.h ../stream_buffer.h | ../bin htslib
	$(CC) $(OPT) $(CFLAGS) -DFASTA=1 -o ../bin/facat seq_cat.c $(LINKING)
../bin/fqcat: seq_cat.c ../seq_file.h ../stream_buffer.h | ../bin htslib
	$(CC) $(OPT) $(CFLAGS) -DFASTQ=1 -o ../bin/fqcat seq_cat.c $(LINKING)
../bin/seqcat: seq_cat.c ../seq_file.h ../stream_buffer.h | ../bin htslib
	$(CC) $(OPT) $(CFLAGS) -DPLAIN=1 -o ../bin/seqcat seq_cat.c $(LINKING)
../bin/revcmp: seq_cat.c ../seq_file.h ../stream_buffer.h | ../bin htslib
	$(CC) $(OPT) $(CFLAGS) -DREVCMP=1 -o ../bin/revcmp seq_cat.c $(LINKING)

../bin/seqrnd: seq_rnd.c
	$(CC) $(OPT) $(CFLAGS) -o ../bin/seqrnd seq_rnd.c

../bin/seqstat: seq_stat.c ../seq_file.h ../stream_buffer.h | ../bin htslib
	$(CC) $(OPT) $(CFLAGS) -o ../bin/seqstat seq_stat.c $(LINKING)

../bin:
	mkdir ../bin

htslib:
	@if [[ '$(HTSLIB)' == '' ]]; \
	then echo "Error: Please pass HTSLIB=... with path to htslib"; exit 1; fi

clean:
	rm -rf ../bin/*  *.dSYM *.greg

.PHONY: all clean htslib
