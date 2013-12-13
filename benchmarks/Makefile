ifdef HTSLIB
	SAMARGS=-I $(HTSLIB)/htslib -L $(HTSLIB) -D_USESAM=1
	SAMLINK=-lhts
endif

CFLAGS=-Wall -Wextra -I .. $(SAMARGS)
LINKING=$(SAMLINK) -lpthread -lz

ifdef DEBUG
	OPT = -O0 --debug -g
else
	OPT = -O3
endif

all: seqtest ktest

seqtest: seq_test.c ../seq_file.h ../stream_buffer.h
	$(CC) $(OPT) $(CFLAGS) -o seqtest seq_test.c $(LINKING)
ktest: ktest.c
	$(CC) $(OPT) -Wall -Wextra -o ktest ktest.c -lz

clean:
	rm -rf seqtest ktest *.dSYM *.greg

.PHONY: all clean htslib
