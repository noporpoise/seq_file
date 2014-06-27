ifdef HTSLIB
	HTSARGS=-I $(HTSLIB)/htslib -D_USESAM=1 $(HTSLIB)/libhts.a -lm
endif

CFLAGS=-Wall -Wextra -std=c99 -pedantic -I ..
LINKING=$(HTSARGS) $(SAMLINK) -lpthread -lz

ifdef DEBUG
	OPT = -O0 --debug -g
else
	OPT = -O3
endif

all: ../bin/facat ../bin/fqcat ../bin/seqcat ../bin/seqstat \
	   ../bin/revcmp ../bin/seqrnd ../bin/dna

../bin/facat: seq_cat.c ../seq_file.h ../stream_buffer.h | ../bin
	$(CC) $(OPT) $(CFLAGS) -DFASTA=1 -o ../bin/facat seq_cat.c $(LINKING)
../bin/fqcat: seq_cat.c ../seq_file.h ../stream_buffer.h | ../bin
	$(CC) $(OPT) $(CFLAGS) -DFASTQ=1 -o ../bin/fqcat seq_cat.c $(LINKING)
../bin/seqcat: seq_cat.c ../seq_file.h ../stream_buffer.h | ../bin
	$(CC) $(OPT) $(CFLAGS) -DPLAIN=1 -o ../bin/seqcat seq_cat.c $(LINKING)
../bin/revcmp: seq_cat.c ../seq_file.h ../stream_buffer.h | ../bin
	$(CC) $(OPT) $(CFLAGS) -DREVCMP=1 -o ../bin/revcmp seq_cat.c $(LINKING)

../bin/seqrnd: seq_rnd.c
	$(CC) $(OPT) $(CFLAGS) -o ../bin/seqrnd seq_rnd.c

../bin/seqstat: seq_stat.c ../seq_file.h ../stream_buffer.h | ../bin
	$(CC) $(OPT) $(CFLAGS) -o ../bin/seqstat seq_stat.c $(LINKING) -lm

../bin/dna: seq_dna.c ../seq_file.h ../stream_buffer.h | ../bin
	$(CC) $(OPT) $(CFLAGS) -o ../bin/dna seq_dna.c $(LINKING) -lm

../bin:
	mkdir ../bin

clean:
	rm -rf ../bin/*  *.dSYM *.greg

.PHONY: all clean
