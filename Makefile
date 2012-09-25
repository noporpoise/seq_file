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

ifndef SAMTOOLS_PATH
	SAMTOOLS_PATH=$(HOME)/bioinf/samtools-0.1.18
endif

# Check mac/linux
UNAME:=$(shell uname)

CFLAGS := $(CFLAGS) -Wall -Wextra -I$(SAMTOOLS_PATH) -I$(STRING_BUF_PATH)

LIB_INCS := -L. -L$(SAMTOOLS_PATH) -L$(STRING_BUF_PATH)
LIB_FLAGS := -lseqfile -lbam -lstrbuf -lz -lm

ifdef ZLIB_PATH
	LIB_INCS := $(LIB_INCS) -L$(ZLIB_PATH)
endif

OBJS = seq_file.o seq_common.o seq_fasta.o seq_fastq.o seq_plain.o seq_sam.o

all: clean $(OBJS)
	ar -csru libseqfile.a $(OBJS)
	$(CC) -o seq_convert $(CFLAGS) $(LIB_INCS) seq_convert.c $(LIB_FLAGS)
	$(CC) -o seq_file_test $(CFLAGS) $(LIB_INCS) seq_file_test.c $(LIB_FLAGS)

clean:
	rm -rf $(OBJS) libseqfile.a seq_convert seq_file_test seq_file.dSYM seq_file.greg

%.o : %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@
