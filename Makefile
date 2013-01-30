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

ifndef HTS_PATH
	HTS_PATH=$(HOME)/bioinf/htslib/
endif

# Argments passed to makefile can only be changed with override
override HTS_PATH := $(HTS_PATH)/htslib/

LIB_STRING_BUF=$(STRING_BUF_PATH)/libstrbuf.a
LIB_HTS=$(HTS_PATH)/libhts.a

CFLAGS := $(CFLAGS) -Wall -Wextra -I $(HTS_PATH) -I $(STRING_BUF_PATH)

LIB_FLAGS := $(LIB_STRING_BUF) $(LIB_HTS) -lpthread -lz -lm

ifdef ZLIB_PATH
	LIB_INCS := $(LIB_INCS) -L $(ZLIB_PATH)
endif

OBJS = seq_file.o seq_common.o seq_fasta.o seq_fastq.o seq_plain.o seq_sam.o

all: clean $(OBJS)
	ar -csru libseqfile.a $(OBJS)
	$(CC) -o seq_convert $(CFLAGS) seq_convert.c libseqfile.a $(LIB_FLAGS)
	$(CC) -o seq_file_test $(CFLAGS) seq_file_test.c libseqfile.a $(LIB_FLAGS)
	cd new_api; make HTS_PATH=$(HTS_PATH)

clean:
	rm -rf $(OBJS) libseqfile.a seq_convert seq_file_test \
	       seq_file.dSYM seq_file.greg seq_convert.dSYM seq_file_test.dSYM
	cd new_api; make clean

.PHONY: all clean

%.o : %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@
