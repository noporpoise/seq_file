seq_file
========

C Library for reading multiple bioinformatics sequence file formats

https://github.com/noporpoise/seq_file

Isaac Turner

turner.isaac@gmail.com

4 July 2012, GPLv3

About
=====

The aim is to provide a C library that allows programs to transparently read
sequence data from multiple file formats, without having to worry about the
format.  Pass a file to seq_file_open() and then read sequences using
seq_next_read() without having to worry about what the file format is.  

Currently supports:
- SAM & BAM
- FASTA (& gzipped fasta)
- FASTQ (& gzipped fastq)
- 'plain' format (one sequence per line [.txt]) (& gzipped plain [.txt.gz])

Also included is a simple example program (seq_convert) that converts between
file formats, for example:

$ seq_convert in.fq.gz out.fa
$ seq_convert in.fa out.fq.gz
$ seq_convert in.bam out.fa.gz
$ seq_convert in.sam out.txt

It can also 'wrap lines' in the output:

$ seq_convert in.fq.gz out.fa 80

Build
=====

Seq_file requires:

sam_tools [http://samtools.sourceforge.net/]
string_buffer [https://github.com/noporpoise/string_buffer]

It also requires zlib, which should already be installed.  

To build the test code and the program seq_convert:

    make STRING_BUF_PATH=path/to/string_buffer/ SAM_PATH=path/to/samtools/

Sometimes the linker can't find your libz.a file (zlib), so you may need to try:

    make STRING_BUF_PATH=path/to/string_buffer/ SAM_PATH=path/to/samtools/ ZLIB_PATH=/dir/with/libz/in/

To call in your own programs, use the following in your Makefile etc.

    LIBS=-lseqfile -lbam -lstrbuf -lz
    INCS=-I$(PATH_TO_SAMTOOLS) -I$(PATH_TO_STRING_BUFFER) -I$(PATH_TO_SEQ_FILE) \
         -L$(PATH_TO_SAMTOOLS) -L$(PATH_TO_STRING_BUFFER) -L$(PATH_TO_SEQ_FILE)
    gcc $(INCS) <your files / args etc.> $(LIBS)

Linking
-------

In order to use in your own

Functions
=========

Opening / Closing
-----------------

Open a file for reading

    SeqFile* seq_file_open(const char* path);

Open a file for writing.  
Returns the number of bytes written or 0 on failure

    SeqFile* seq_file_open_write(const char* file_path, const SeqFileType file_type,
                                 const char gzip, const unsigned long line_wrap);

Close a file. If open for writing: writes newline to file and returns 1 on
success, otherwise returns 0.

    size_t seq_file_close(SeqFile *sf);

Reading
-------

Get the next sequence read.  Returns 1 on success 0 if no more to read

    char seq_next_read(SeqFile *sf);

Get the name of the next read

    const char* seq_get_read_name(SeqFile *sf);

Get this read index -- starts from 0

    unsigned long seq_get_read_index(SeqFile *sf);

Get the distance into this read that we have read

    unsigned long seq_get_base_offset(SeqFile *sf);

Get the distance into this read's quality scores that we have read

    unsigned long seq_get_qual_offset(SeqFile *sf);

If seq_next_read() returned 1 and seq_read_base() is now returning 0,
seq_get_length() will now report the correct read length

    unsigned long seq_get_length(SeqFile *sf);

Read a single base from the current read
Returns 1 on success, 0 if no more quality scores or run out of bases

    char seq_read_base(SeqFile *sf, char *c);

Read a single quality score from the current read
Returns 1 on success, 0 if no more quality scores or run out of bases

    char seq_read_qual(SeqFile *sf, char *c);

Read `k` bases of the current read.  `str` must be at least k+1 bytes long.
Returns 1 on success, 0 otherwise

    char seq_read_k_bases(SeqFile *sf, char* str, int k);
    char seq_read_k_quals(SeqFile *sf, char* str, int k);

Read all remaining bases from the current read.
Returns 1 on success, 0 otherwise.

    char seq_read_all_bases(SeqFile *sf, StrBuf *sbuf);
    char seq_read_all_quals(SeqFile *sf, StrBuf *sbuf);

Writing
-------

Returns 1 if open for writing, 0 otherwise

    char seq_is_open_for_write(const SeqFile *sf);

Write various values to the file.  Must be called in order of: name, seq, qual.
Any of `name`, `seq`, `quals` may be NULL.  (These API calls may change).  

    size_t seq_file_write_name(SeqFile *sf, const char *name);
    size_t seq_file_write_seq(SeqFile *sf, const char *seq);
    size_t seq_file_write_qual(SeqFile *sf, const char *qual);

Get file details
----------------

Guess a filetype from path

    void seq_guess_filetype_from_path(const char* path,
                                      SeqFileType* file_type,
                                      char* zipped);

Various methods for getting file type

    SeqFileType seq_file_get_type(const SeqFile *sf);
    const char* seq_file_get_type_str(const SeqFile  *sf);
    const char* seq_file_type_str(const SeqFileType file_type, const char zipped);

Get file path

    const char* seq_get_path(const SeqFile *sf);

Get the number of bases read/written so far

    unsigned long seq_num_bases_passed(const SeqFile *sf);

Whether or not this file has quality scores.  Returns 1 if file has quality
scores, 0 otherwise.  Note: quality scores may still all be set to 'null' e.g.
??? or ### etc.

    char seq_has_quality_scores(const SeqFile *sf);

Development
===========

Please get in touch if you find bugs, have questions or can suggest features.  
Isaac <turner.isaac@gmail.com>

To run some basic test:

    for f in test/*; do echo $f >&2; ./seq_file_test $f; done > test/tests.out

TODO:
seq_file
 * Retrieve mate pair data etc from a sam/bam
 * Add support for sra
 
seq_convert
 * Add pair-end awareness: bam to 2 fastq files


License
=======

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
