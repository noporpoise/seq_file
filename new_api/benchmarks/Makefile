ifndef HTSLIB
	HTSLIB=$(HOME)/bioinf/htslib/
endif
# Argments passed to makefile can only be changed with override
override HTSLIB:=$(HTSLIB)/htslib/

CFLAGS=-Wall -Wextra -I .. -I $(HTSLIB) -L $(HTSLIB)
LINKING=-lhts -lpthread -lz

all: seqtest ktest

seqtest: seq_test.c ../seq_file.h ../buffered_input.h
	$(CC) -O3 $(CFLAGS) -o seqtest seq_test.c $(LINKING)
ktest: ktest.c kseq.h
	$(CC) -O3 -o ktest ktest.c -lz

clean:
	rm -rf seqtest ktest

.PHONY: all clean
