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
SequenceFile* seq_file_open(const char* path);
SequenceFile* seq_file_open_filetype(const char* file_path,
                                     const SequenceFileType file_type);

void seq_file_close(SequenceFile* file);

SequenceFileType seq_file_get_type(const SequenceFile* file);
const char* seq_file_get_type_str(const SequenceFile* file);

const char* seq_file_get_path(const SequenceFile* file);

// Returns 0 if at end of file; 1 otherwise
char seq_file_read(SequenceFile* file, Sequence* sequence);

//
// Sequence Kmer Reader
//
SequenceKmerReader* seq_file_kmer_reader_open(const char* path,
                                              const int kmer_size);

void seq_file_kmer_reader_close(SequenceKmerReader* reader);

// Returns 0 if at end of file; 1 otherwise
char seq_file_read_kmer(SequenceKmerReader* reader, Sequence* sequence);

#endif
