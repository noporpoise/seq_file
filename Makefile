ifndef CC
  CC = gcc
endif

ifdef DEBUG
	CFLAGS := -DDEBUG=1 --debug -g
else
	CFLAGS := -O3
endif

ifndef STRING_BUF_PATH
	STRING_BUF_PATH=$(HOME)/c/libs/string_buffer
endif

ifndef SAM_PATH
	SAMPATH=$(HOME)/bioinf/samtools-0.1.18
endif

# Check mac/linux
UNAME:=$(shell uname)

CFLAGS := $(CFLAGS) -Wall -Wextra -I$(SAMPATH) -I$(STRING_BUF_PATH)

ifeq ($(CC),gcc)
  ifeq ($(UNAME), Darwin)
    CFLAGS := $(CFLAGS) -fnested-functions
  endif
endif

LIB_INCS = -L$(SAMPATH) -L$(STRING_BUF_PATH)
LIB_FLAGS := -lbam -lm -lz -lstrbuf

all:
	$(CC) $(CFLAGS) -o seq_file.o -c seq_file.c
	ar -csru libseqfile.a seq_file.o
	$(CC) -o seq_convert $(CFLAGS) $(LIB_INCS) seq_convert.c seq_file.o $(LIB_FLAGS)
	$(CC) -o seq_file_test $(CFLAGS) $(LIB_INCS) seq_file_test.c seq_file.o $(LIB_FLAGS)

clean:
	if test -e seq_file.o; then rm seq_file.o; fi
	if test -e libseqfile.a; then rm libseqfile.a; fi
	if test -e seq_convert; then rm seq_convert; fi
	if test -e seq_file_test; then rm seq_file_test; fi
	for file in $(wildcard *.dSYM); do rm -r $$file; done
	for file in $(wildcard *.greg); do rm $$file; done
