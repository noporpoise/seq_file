ifndef CC
	CC=gcc
endif

HTS=$(HTSLIB)/htslib/
CFLAGS=-Wall -Wextra

ifdef DEBUG
	OPT = -O0 --debug -g
else
	OPT = -O3
endif

all: seqtest ktest

seqtest: seq_test.c ../seq_file.h ../buffered_input.h | htslib
	$(CC) $(OPT) $(CFLAGS) -I.. -I$(HTS) -L$(HTS) -o seqtest seq_test.c -lhts -lpthread -lz
ktest: ktest.c | htslib
	$(CC) $(OPT) $(CFLAGS) -I$(HTS) -o ktest ktest.c -lz

htslib:
	if [[ '$(HTSLIB)' == '' ]]; \
	then echo "Error: Please pass HTSLIB=... with path to htslib"; exit 1; fi

clean:
	rm -rf seqtest ktest *.dSYM *.greg

.PHONY: all clean htslib
