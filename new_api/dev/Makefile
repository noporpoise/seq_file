ifndef CC
	CC=gcc
endif

CFLAGS=-Wall -Wextra -I .. -I $(HTSLIB)/htslib/ -L $(HTSLIB)/htslib/
LINKING=-lhts -lpthread -lz

ifdef DEBUG
	OPT = -O0 --debug -g
else
	OPT = -O3
endif

all: buftest seqcmdline

buftest: buf_test.c ../buffered_input.h | htslib
	$(CC) $(OPT) -I .. -o buftest buf_test.c -lz 
seqcmdline: seq_cmdline.c ../seq_file.h ../buffered_input.h | htslib
	$(CC) $(OPT) $(CFLAGS) -o seqcmdline seq_cmdline.c $(LINKING)

htslib:
	@if [[ '$(HTSLIB)' == '' ]]; \
	then echo "Error: Please pass HTSLIB=... with path to htslib"; exit 1; fi

clean:
	rm -rf buftest seqcmdline *.dSYM *.greg

.PHONY: all clean htslib
