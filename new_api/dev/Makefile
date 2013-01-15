ifndef HTSLIB
	HTSLIB=$(HOME)/bioinf/htslib/
endif
# Argments passed to makefile can only be changed with override
override HTSLIB:=$(HTSLIB)/htslib/

CFLAGS=-I .. -I $(HTSLIB) -L $(HTSLIB)
LINKING=-lhts -lpthread -lz

all: btest seqcmdline seqstat

btest: buffered_input.c ../buffered_input.h
	$(CC) -O3 -I .. -o btest buffered_input.c -lz
seqcmdline: seq_cmdline.c ../seq_file.h ../buffered_input.h
	$(CC) -O3 $(CFLAGS) -o seqcmdline seq_cmdline.c $(LINKING)
seqstat: seq_stat.c ../seq_file.h ../buffered_input.h
	$(CC) $(CFLAGS) -o seqstat seq_stat.c $(LINKING)

clean:
	rm -rf btest seqcmdline seqstat

.PHONY: all clean
