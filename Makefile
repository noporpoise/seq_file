ifdef DEBUG
	FLAGS=-DDEBUG=1 --debug
endif

LIB_PATH=$(HOME)/c/libs/

STRING_BUF_PATH=$(LIB_PATH)/string_buffer
SAMPATH=$(HOME)/bioinf/samtools-0.1.18

# Check mac/linux
UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
	FLAGS := $(FLAGS) -fnested-functions
endif

FLAGS := $(FLAGS) -Wall -I$(SAMPATH) -L$(SAMPATH) -lbam -lm -lz

all:
	gcc -o seq_reader_test $(FLAGS) -I $(STRING_BUF_PATH) \
	seq_reader_test.c seq_reader.c \
	$(STRING_BUF_PATH)/string_buffer.c

clean:
	if test -e seq_reader_test; then rm seq_reader_test; fi
	for file in $(wildcard *.dSYM); do rm -r $$file; done
	for file in $(wildcard *.greg); do rm $$file; done
