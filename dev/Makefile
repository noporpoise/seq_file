ifdef HTSLIB
	SAMARGS=-I $(HTSLIB)/htslib -L $(HTSLIB) -D_USESAM=1
	SAMLINK=-lhts -lm
endif

CFLAGS=-Wall -Wextra -I .. $(SAMARGS)
LINKING=$(SAMLINK) -lpthread -lz

ifdef DEBUG
	OPT = -O0 --debug -g
else
	OPT = -O3
endif

all: seqcmdline

seqcmdline: seq_cmdline.c ../seq_file.h ../stream_buffer.h
	$(CC) $(OPT) $(CFLAGS) -o seqcmdline seq_cmdline.c $(LINKING)

clean:
	rm -rf seqcmdline *.dSYM *.greg

.PHONY: all clean
