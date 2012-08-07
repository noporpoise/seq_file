ifdef DEBUG
	CFLAGS := -DDEBUG=1 --debug
else
	CFLAGS := -O3
endif

LIB_PATH=$(HOME)/c/libs

STRING_BUF_PATH=$(LIB_PATH)/string_buffer
SAMPATH=$(HOME)/bioinf/samtools-0.1.18

# Check mac/linux
UNAME:=$(shell uname)

ifeq ($(UNAME), Darwin)
	CFLAGS := $(CFLAGS) -fnested-functions
endif

CFLAGS := $(CFLAGS) -Wall -Wextra -L$(SAMPATH) -I$(SAMPATH) \
		  -I$(STRING_BUF_PATH) -L$(STRING_BUF_PATH)

LIB_FLAGS := -lbam -lm -lz -lstrbuf

all:
	gcc $(CFLAGS) -o seq_file.o -c seq_file.c
	ar -csru libseqfile.a seq_file.o
	gcc -o seq_convert $(CFLAGS) seq_convert.c seq_file.o $(LIB_FLAGS)
	gcc -o seq_file_test $(CFLAGS) seq_file_test.c seq_file.o $(LIB_FLAGS)

clean:
	if test -e seq_file.o; then rm seq_file.o; fi
	if test -e libseqfile.a; then rm libseqfile.a; fi
	if test -e seq_convert; then rm seq_convert; fi
	if test -e seq_file_test; then rm seq_file_test; fi
	for file in $(wildcard *.dSYM); do rm -r $$file; done
	for file in $(wildcard *.greg); do rm $$file; done
