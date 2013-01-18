ifndef HTSLIB
	HTSLIB=$(HOME)/bioinf/htslib/
endif
# Argments passed to makefile can only be changed with override
override HTSLIB:=$(HTSLIB)/htslib/

CFLAGS=-O3 -I .. -I $(HTSLIB) -L $(HTSLIB)
LINKING=-lhts -lpthread -lz

all: ../bin/facat ../bin/fqcat ../bin/seqcat

../bin/facat: ../bin seq_cat.c ../seq_file.h ../buffered_input.h
	$(CC) -DFASTA=1 $(CFLAGS) -o ../bin/facat seq_cat.c $(LINKING)
../bin/fqcat: ../bin seq_cat.c ../seq_file.h ../buffered_input.h
	$(CC) -DFASTQ=1 $(CFLAGS) -o ../bin/fqcat seq_cat.c $(LINKING)
../bin/seqcat: ../bin seq_cat.c ../seq_file.h ../buffered_input.h
	$(CC) -DPLAIN=1 $(CFLAGS) -o ../bin/seqcat seq_cat.c $(LINKING)

../bin:
	mkdir -p ../bin

clean:
	rm -rf ../bin/*

.PHONY: all clean
