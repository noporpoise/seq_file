/*
 seq_reader.c
 project: seq_reader
 author: Isaac Turner <turner.isaac@gmail.com>
 Copyright (C) 20-June-2012
 
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
 */

#ifndef SEQ_HEADER_SEEN
#define SEQ_HEADER_SEEN

#include "string_buffer.h"

typedef struct SequenceFile SequenceFile;
typedef struct SequenceKmerReader SequenceKmerReader;
typedef enum SequenceFileType SequenceFileType;

typedef struct
{
  STRING_BUFFER *name, *seq, *qual;
  int offset;
  char valid_qual;
} Sequence;

enum SequenceFileType
{
  SEQ_UNKNOWN = 0, SEQ_FASTA = 1, SEQ_FASTQ = 2, SEQ_PLAIN = 3,
  SEQ_SAM = 4, SEQ_BAM = 5,
};

//
// Sequence struct
//

// For creating/destroying struct for result
Sequence* seq_create();
void seq_reset(Sequence* sequence);
void seq_free(Sequence* sequence);

//
// Sequence File reader
//
// Open to read
SequenceFile* seq_file_open(const char* file_path);

// Open to read and force a particular file type
SequenceFile* seq_file_open_filetype(const char* file_path,
                                     const SequenceFileType file_type);

// Open to write
// file_type must be FASTA, FASTQ
// set gzip to != 0 to turn on gzipping output
// if line_wrap is != 0, sequence lines are wrapped
SequenceFile* seq_file_open_write(const char* file_path,
                                  const SequenceFileType file_type,
                                  const char gzip,
                                  const int line_wrap);

// Close
void seq_file_close(SequenceFile* file);

// Guess a filetype from path
void seq_file_guess_filetype_from_path(const char* path,
                                       SequenceFileType* file_type,
                                       char* zipped);

// Get file type
SequenceFileType seq_file_get_type(const SequenceFile* file);
const char* seq_file_get_type_str(const SequenceFile* file);
const char* seq_file_convert_type_str(const SequenceFileType file_type,
                                      const char zipped);

// Get path
const char* seq_file_get_path(const SequenceFile* file);

// Read a whole sequence from the file
// Returns 1 on success, 0 otherwise
char seq_file_read(SequenceFile* file, Sequence* sequence);

// Get the current entry number
// 0 when nothing read/written
// 1 after first read/write etc.
long seq_file_entry_number(SequenceFile *file);

// Write a sequence
// returns number of bytes written (0 on failure)
unsigned long seq_file_write(SequenceFile *file, Sequence *seq);

//
// Sequence Kmer Reader
//
// Open and close:
SequenceKmerReader* kmer_reader_open(const char* path, const int kmer_size);
void kmer_reader_close(SequenceKmerReader* reader);

// To read:

// Returns 0 if at end of file; 1 otherwise
char kmer_reader_read_kmer(SequenceKmerReader* reader, Sequence* sequence);

// Return next char or -1 if end of read entry
int kmer_reader_read_char(SequenceKmerReader* reader);

long kmer_reader_entry_number(SequenceKmerReader* reader);

#endif
