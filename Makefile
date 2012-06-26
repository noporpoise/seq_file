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

CFLAGS := $(CFLAGS) -Wall -I$(SAMPATH) -L$(SAMPATH) -I $(STRING_BUF_PATH)
LIB_FLAGS := -lbam -lm -lz

all:
	gcc -o seq_reader.o $(CFLAGS) -c seq_reader.c
	gcc -o seq_convert $(CFLAGS) $(LIB_FLAGS) \
	seq_convert.c seq_reader.o $(STRING_BUF_PATH)/string_buffer.c
	gcc -o seq_reader_test $(CFLAGS) $(LIB_FLAGS) \
	seq_reader_test.c seq_reader.o $(STRING_BUF_PATH)/string_buffer.c

clean:
	if test -e seq_reader.o; then rm seq_reader.o; fi
	if test -e seq_convert; then rm seq_convert; fi
	if test -e seq_reader_test; then rm seq_reader_test; fi
	for file in $(wildcard *.dSYM); do rm -r $$file; done
	for file in $(wildcard *.greg); do rm $$file; done
