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

CFLAGS := $(CFLAGS) -Wall -Wextra
CFLAGS := $(CFLAGS) -I$(SAMPATH) -L$(SAMPATH) -I $(STRING_BUF_PATH)

LIB_FLAGS := -lbam -lm -lz

all:
	gcc -o seq_file.o $(CFLAGS) -c seq_file.c
	gcc -o seq_convert $(CFLAGS) $(LIB_FLAGS) \
	seq_convert.c seq_file.o $(STRING_BUF_PATH)/string_buffer.c
	gcc -o seq_file_test $(CFLAGS) $(LIB_FLAGS) \
	seq_file_test.c seq_file.o $(STRING_BUF_PATH)/string_buffer.c

clean:
	if test -e seq_file.o; then rm seq_file.o; fi
	if test -e seq_convert; then rm seq_convert; fi
	if test -e seq_file_test; then rm seq_file_test; fi
	for file in $(wildcard *.dSYM); do rm -r $$file; done
	for file in $(wildcard *.greg); do rm $$file; done
