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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "seq_reader.h"
#include "sam.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))

#define is_base_char(x) ((x) == 'a' || (x) == 'A' || \
                         (x) == 'c' || (x) == 'C' || \
                         (x) == 'g' || (x) == 'G' || \
                         (x) == 't' || (x) == 'T' || \
                         (x) == 'n' || (x) == 'N')

int8_t seq_comp_table[16] = {0,8,4,12,2,10,9,14,1,6,5,13,3,11,7,15};

const char* seq_file_types[6]
  = {"Unknown", "FASTA", "FASTQ", "Plain", "SAM", "BAM"};

const char* seq_file_types_zipped[6]
  = {"Unknown(zipped)", "FASTA(zipped)", "FASTQ(zipped)", "Plain(zipped)",
     "SAM", "BAM"};

struct SequenceFile
{
  const char *path;

  //union {
  gzFile *gz_file; // for reading
  FILE *plain_file; // for writing
  samfile_t *sam_file; // For reading sam/bams
  //} f;

  // For reading sam/bams
  bam1_t *bam;

  // have we seen a '>' at the start of a line in a fasta file?
  char read_line_start;

  unsigned long entry_number; // 1-based, 0 means haven't started reading

  enum SequenceFileType file_type;
  
  // Writing preferences
  int line_wrap;
};

struct SequenceKmerReader
{
  SequenceFile* file;
  int kmer_size;

  // points to the character just returned
  unsigned long offset;

  // need to read in whole entry for sam/bam and fastq
  Sequence* whole_entry;
};


//
// Sequence struct
//

// For creating/destroying struct for result
Sequence* seq_create()
{
  Sequence* sequence = (Sequence*) malloc(sizeof(Sequence));

  sequence->name = string_buff_init(100);
  sequence->seq = string_buff_init(100);
  sequence->qual = string_buff_init(100);

  sequence->offset = 0;
  sequence->valid_qual = 0;

  return sequence;
}

void seq_reset(Sequence* sequence)
{
  string_buff_reset(sequence->name);
  string_buff_reset(sequence->seq);
  string_buff_reset(sequence->qual);
  sequence->offset = 0;
  sequence->valid_qual = 0;
}

void seq_free(Sequence* sequence)
{
  string_buff_free(sequence->name);
  string_buff_free(sequence->seq);
  string_buff_free(sequence->qual);

  free(sequence);
}


//
// Sequence File reader
//

void seq_file_guess_filetype_from_path(const char* path,
                                       SequenceFileType* file_type,
                                       char* zipped)
{
  size_t path_len = strlen(path);

  if(strcasecmp(".sam", path+path_len - 4) == 0)
  {
    *file_type = SEQ_SAM;
  }
  else if(strcasecmp(".bam", path+path_len - 4) == 0)
  {
    *file_type = SEQ_BAM;
  }
  else if(strcasecmp(".fa", path+path_len - 3) == 0 ||
          strcasecmp(".fasta", path+path_len - 6) == 0)
  {
    *zipped = 0;
    *file_type = SEQ_FASTA;
  }
  else if(strcasecmp(".fa.gz", path+path_len - 6) == 0 ||
          strcasecmp(".fasta.gz", path+path_len - 9) == 0)
  {
    *zipped = 1;
    *file_type = SEQ_FASTA;
  }
  else if(strcasecmp(".fq", path+path_len - 3) == 0 ||
          strcasecmp(".fastq", path+path_len - 6) == 0)
  {
    *zipped = 0;
    *file_type = SEQ_FASTQ;
  }
  else if(strcasecmp(".fq.gz", path+path_len - 6) == 0 ||
          strcasecmp(".fastq.gz", path+path_len - 9) == 0)
  {
    *zipped = 1;
    *file_type = SEQ_FASTQ;
  }
  else
  {
    *file_type = SEQ_UNKNOWN;
  }
}

void _set_seq_filetype(SequenceFile* file)
{
  size_t path_len = strlen(file->path);

  if(path_len > 4)
  {
    if(strcasecmp(".sam", file->path+path_len - 4) == 0)
    {
      file->sam_file = samopen(file->path, "r", 0);
      file->bam = bam_init1();
      file->file_type = SEQ_SAM;
      return;
    }
    else if(strcasecmp(".bam", file->path+path_len - 4) == 0)
    {
      file->sam_file = samopen(file->path, "rb", 0);
      file->bam = bam_init1();
      file->file_type = SEQ_BAM;
      return;
    }
  }

  // Open file for the first time
  if(strcmp(file->path, "-") == 0)
  {
    file->gz_file = gzdopen(fileno(stdin), "r");
  }
  else
  {
    file->gz_file = gzopen(file->path, "r");
  }

  int first_char;

  do
  {
    first_char = gzgetc(file->gz_file);
  } while (first_char != -1 && (first_char == '\n' || first_char == '\r'));

  if(first_char == -1)
  {
    fprintf(stderr, "seq_reader.c warning: empty sequence file\n");
    return;
  }
  else if(first_char == '>')
  {
    // Reading FASTA
    file->file_type = SEQ_FASTA;
    file->read_line_start = 1;
  }
  else if(first_char == '@')
  {
    // Reading FASTQ
    file->file_type = SEQ_FASTQ;
    file->read_line_start = 1;
  }
  else if(is_base_char(first_char))
  {
    file->file_type = SEQ_PLAIN;
    file->read_line_start = first_char;
  }
  else
  {
    fprintf(stderr, "seq_reader.c warning: unknown filetype starting '%c'\n",
            first_char);
  }
}

SequenceFile* seq_file_open(const char* file_path)
{
  SequenceFile* file = (SequenceFile*) malloc(sizeof(SequenceFile));
  file->path = file_path;
  file->gz_file = NULL;
  file->file_type = SEQ_UNKNOWN;
  file->sam_file = NULL;
  file->bam = NULL;
  file->line_wrap = 0;
  file->read_line_start = 0;
  file->entry_number = 0;

  _set_seq_filetype(file);

  if(file->file_type == SEQ_UNKNOWN)
  {
    free(file);
    file = NULL;
  }

  return file;
}

SequenceFile* seq_file_open_filetype(const char* file_path,
                                     const SequenceFileType file_type)
{
  SequenceFile* file = (SequenceFile*) malloc(sizeof(SequenceFile));
  file->path = file_path;
  file->gz_file = NULL;
  file->file_type = file_type;
  file->sam_file = NULL;
  file->bam = NULL;
  file->line_wrap = 0;
  file->read_line_start = 0;
  file->entry_number = 0;

  switch(file_type)
  {
    case SEQ_SAM:
      file->sam_file = samopen(file->path, "r", 0);
      file->bam = bam_init1();
      break;
    case SEQ_BAM:
      file->sam_file = samopen(file->path, "rb", 0);
      file->bam = bam_init1();
      break;
    case SEQ_FASTA:
    case SEQ_FASTQ:
    case SEQ_PLAIN:
      if(strcmp(file->path, "-") == 0)
      {
        file->gz_file = gzdopen(fileno(stdin), "r");
      }
      else
      {
        file->gz_file = gzopen(file->path, "r");
      }
      break;
    default:
      fprintf(stderr, "seq_reader.c warning: invalid SequenceFileType in "
                      "function seq_file_open_filetype()\n");
      free(file);
      return NULL;
  }

  return file;
}

// Open to write
// file_type must be FASTA, FASTQ
// set gzip to != 0 to turn on gzipping output
// if line_wrap is != 0, sequence lines are wrapped
SequenceFile* seq_file_open_write(const char* file_path,
                                  const SequenceFileType file_type,
                                  const char gzip,
                                  const int line_wrap)
{
  if(file_type == SEQ_SAM || file_type == SEQ_BAM)
  {
    fprintf(stderr, "seq_reader.c error: cannot write to a SAM or BAM file\n");
    return NULL;
  }
  else if(file_type == SEQ_UNKNOWN)
  {
    fprintf(stderr, "seq_reader.c error: cannot open file type SEQ_UNKNOWN\n");
    return NULL;
  }
  else if(file_type == SEQ_PLAIN && line_wrap != 0)
  {
    fprintf(stderr, "seq_reader.c warning: cannot set line wrap with 'plain' "
                    "sequence format\n");
  }

  gzFile* gz_file = NULL;
  FILE* plain_file = NULL;

  if(gzip)
  {
    if((gz_file = gzopen(file_path, "w")) == NULL)
    {
      return NULL;
    }
  }
  else
  {
    if((plain_file = fopen(file_path, "w")) == NULL)
    {
      return NULL;
    }
  }

  SequenceFile* file = (SequenceFile*) malloc(sizeof(SequenceFile));
  file->path = file_path;
  file->gz_file = gz_file;
  file->plain_file = plain_file;
  file->file_type = file_type;
  file->sam_file = NULL;
  file->bam = NULL;
  file->line_wrap = (file_type == SEQ_PLAIN && line_wrap != 0) ? 0 : line_wrap;
  file->read_line_start = 0;
  file->entry_number = 0;

  return file;
}


// Close
void seq_file_close(SequenceFile* file)
{
  if(file->gz_file != NULL)
  {
    gzclose(file->gz_file);
  }

  if(file->file_type == SEQ_SAM || file->file_type == SEQ_BAM)
  {
    bam_destroy1(file->bam);
    samclose(file->sam_file);
  }

  free(file);
}

SequenceFileType seq_file_get_type(const SequenceFile* file)
{
  return file->file_type;
}

const char* seq_file_get_type_str(const SequenceFile* file)
{
  return seq_file_types[file->file_type];
}

const char* seq_file_convert_type_str(const SequenceFileType file_type,
                                      const char zipped)
{
  return zipped ? seq_file_types_zipped[file_type] : seq_file_types[file_type];
}

// Get a pointer to the file path
const char* seq_file_get_path(const SequenceFile* file)
{
  return file->path;
}

// Read an entry from a FASTQ file
// Returns 1 if a header is read, 0 otherwise
// Prints errors to stderr if there are syntax issues
char _read_fastq_entry(SequenceFile* file, Sequence* sequence)
{
  // Read until we have a header line
  int c;

  if(!file->read_line_start)
  {
    do
    {
      c = gzgetc(file->gz_file);
    } while(c != -1 && (c == '\n' || c == '\r'));

    if(c == -1)
    {
      return 0;
    }
    else if(c != '@')
    {
      fprintf(stderr, "seq_reader.c warning: FASTQ header does not begin "
                      "with '@'\n");
      return 0;
    }
  }

  file->read_line_start = 0;

  int first_char;

  // Read name
  if(string_buff_gzreadline(sequence->name, file->gz_file) == 0)
  {
    fprintf(stderr, "seq_reader.c warning: FASTQ missing name line "
                    "[read: '%s']\n", sequence->name->buff);
    return 0;
  }

  string_buff_chomp(sequence->name);

  // Read sequence
  while((first_char = gzgetc(file->gz_file)) != -1 && first_char != '+')
  {
    if(first_char != '\n' && first_char != '\r')
    {
      string_buff_append_char(sequence->seq, first_char);
      string_buff_gzreadline(sequence->seq, file->gz_file);
      string_buff_chomp(sequence->seq);
    }
  }

  if(string_buff_strlen(sequence->seq) == 0)
  {
    fprintf(stderr, "seq_reader.c warning: FASTQ missing sequence line "
                    "[read: '%s']\n", sequence->name->buff);
    return 0;
  }
  else if(first_char == -1)
  {
    fprintf(stderr, "seq_reader.c warning: FASTQ ended prematurely "
                    "[read: '%s']\n", sequence->name->buff);
    return 0;
  }
  else if(first_char != '+')
  {
    fprintf(stderr, "seq_reader.c warning: FASTQ separator '+' missing "
                    "[read: '%s']\n", sequence->name->buff);
    return 0;
  }

  // Read quality line
  while((first_char = gzgetc(file->gz_file)) != -1 &&
        string_buff_strlen(sequence->qual) < string_buff_strlen(sequence->seq))
  {
    if(first_char != '\n' && first_char != '\r')
    {
      string_buff_append_char(sequence->qual, first_char);
      string_buff_gzreadline(sequence->qual, file->gz_file);
      string_buff_chomp(sequence->qual);
    }
  }
  
  if(first_char == '@')
  {
    file->read_line_start = 1;
  }

  if(string_buff_strlen(sequence->qual) == 0)
  {
    fprintf(stderr, "seq_reader.c warning: FASTQ is missing a quality line "
                    "[read: '%s']\n", sequence->name->buff);
  }
  else if(string_buff_strlen(sequence->qual) !=
          string_buff_strlen(sequence->seq))
  {
    fprintf(stderr, "seq_reader.c warning: FASTQ sequence and quality lines "
                    "are not the same length [read: '%s']\n",
                    sequence->name->buff);
    //fprintf(stderr, "Sequence: '%s'\n", sequence->seq->buff);
    //fprintf(stderr, "Quality : '%s'\n", sequence->qual->buff);
  }
  else
  {
    sequence->valid_qual = 1;
  }

  return 1;
}

// Read an entry from a FASTA file
// Returns 1 if a header is read, 0 otherwise
char _read_fasta_entry(SequenceFile* file, Sequence* sequence)
{
  // Read until we have a header line
  int c;

  if(!file->read_line_start)
  {
    do
    {
      c = gzgetc(file->gz_file);
    } while(c != -1 && (c == '\n' || c == '\r'));

    if(c == -1)
    {
      return 0;
    }
    else if(c != '>')
    {
      fprintf(stderr, "seq_reader.c warning: FASTA header does not begin "
                      "with '>' [read name: '%s']\n", sequence->name->buff);
    }
  }

  file->read_line_start = 0;

  string_buff_gzreadline(sequence->name, file->gz_file);
  string_buff_chomp(sequence->name);

  // Check line doesn't begin with '>'
  int first_char;

  while((first_char = gzgetc(file->gz_file)) != -1)
  {
    if(first_char == '>')
    {
      // Done
      file->read_line_start = 1;
      break;
    }
    else if(first_char != '\n' && first_char != '\r')
    {
      // Push char onto string
      string_buff_append_char(sequence->seq, first_char);
      string_buff_chomp(sequence->seq);

      // Read the rest of the line
      if(string_buff_gzreadline(sequence->seq, file->gz_file) == 0)
      {
        break;
      }

      string_buff_chomp(sequence->seq);
    }
  }

  sequence->valid_qual = 0;

  return 1;
}

// Read an entry from a SAM/BAM file
// Returns 1 if an entry is read, 0 otherwise
char _read_bam_entry(SequenceFile* file, Sequence* sequence)
{
  if(samread(file->sam_file, file->bam) >= 0)
  {
    uint8_t *seq; // samtools gives us seq pointer
    char* str;
    int i;

    int qlen = file->bam->core.l_qseq;

    // Get reverse
    char is_reversed = file->bam->core.flag & 16;

    //
    // Get name
    //
    string_buff_append_str(sequence->name, bam1_qname(file->bam));

    //
    // Get sequence
    //
    string_buff_ensure_capacity(sequence->seq, qlen);
    sequence->seq->len = qlen;
    str = sequence->seq->buff;
    str[qlen] = '\0';

    seq = bam1_seq(file->bam);

    // read in and complement (if needed)
    for(i = 0; i < qlen; i++)
    {
      int8_t b = bam1_seqi(seq, i);
      str[i] = bam_nt16_rev_table[is_reversed ? seq_comp_table[b] : b];
    }

    // reverse
    if(is_reversed)
    {
      string_buff_reverse(sequence->seq);
    }

    //
    // Get quality
    //
    string_buff_ensure_capacity(sequence->qual, qlen);
    sequence->qual->len = qlen;
    str = sequence->qual->buff;
    str[qlen] = '\0';

    seq = bam1_qual(file->bam);

    if(is_reversed)
    {
      // Read-in backwards to do reverse
      for (i = 0; i < qlen; i++)
        str[qlen - i - 1] = 33 + seq[i];
    }
    else
    {
      // Read-in forward
      for (i = 0; i < qlen; i++)
        str[i] = 33 + seq[i];
    }

    sequence->valid_qual = 0;

    // Check if the quality string has any values set    
    for (i = 0; i < qlen; i++)
    {
      if(str[i] != '?')
      {
        sequence->valid_qual = 1;
        break;
      }
    }

    return 1;
  }

  return 0;
}

// Read an entry from a SAM/BAM file
// Returns 1 if a line is read, 0 otherwise
char _read_plain_entry(SequenceFile* file, Sequence* sequence)
{
  if(file->read_line_start)
  {
    string_buff_append_char(sequence->seq, file->read_line_start);
    file->read_line_start = 0;
  }

  // If we don't read any bytes, return 0
  if(string_buff_gzreadline(sequence->seq, file->gz_file) == 0)
  {
    return 0;
  }

  string_buff_chomp(sequence->seq);
  sequence->valid_qual = 0;

  return 1;
}

// Returns 0 if at end of file; 1 otherwise
char seq_file_read(SequenceFile* file, Sequence* sequence)
{
  seq_reset(sequence);

  if(file->gz_file != NULL && gzeof(file->gz_file))
  {
    return 0;
  }

  char success;

  switch(file->file_type)
  {
    case SEQ_SAM:
    case SEQ_BAM:
      success = _read_bam_entry(file, sequence);
      break;
    case SEQ_FASTA:
      success = _read_fasta_entry(file, sequence);
      break;
    case SEQ_FASTQ:
      success = _read_fastq_entry(file, sequence);
      break;
    case SEQ_PLAIN:
      success = _read_plain_entry(file, sequence);
      break;
    default:
      fprintf(stderr, "seq_reader.c warning: "
                      "tried to read from unknown filetype\n");
      return 0;
  }

  if(success)
  {
    file->entry_number++;
  }

  sequence->offset = 0;

  return success;
}

//
// Write
//

inline unsigned long _write(SequenceFile* file, const char* str, const unsigned long len)
{
  if(file->plain_file != NULL)
  {
    return fwrite(str, sizeof(char), len, file->plain_file);
  }
  else if(gzwrite(file->gz_file, str, len))
  {
    return len;
  }

  return 0;
}

inline unsigned long _write_wrapped(STRING_BUFFER* sb, SequenceFile *file)
{
  unsigned long len = 0, offset, line_len;

  for(offset = 0; offset < sb->len; offset += file->line_wrap)
  {
    line_len = MIN(sb->len - offset, file->line_wrap);

    if(!_write(file, sb->buff + offset, line_len) || !_write(file, "\n", 1))
    {
      return 0;
    }

    len += line_len + 1;
  }

  if(sb->len == 0)
  {
    if(!_write(file, "\n", 1))
    {
      return 0;
    }

    len += 1;
  }

  return len;
}

unsigned long _write_fasta_entry(SequenceFile *file, Sequence *seq)
{
  unsigned long len = 0;

  if(!_write(file, ">", 1) ||
     (seq->name->len != 0 && !_write(file, seq->name->buff, seq->name->len)) ||
     !_write(file, "\n", 1))
  {
    return 0;
  }

  len += 1 + seq->name->len + 1;

  if(file->line_wrap == 0)
  {
    if((seq->seq->len != 0 && !_write(file, seq->seq->buff, seq->seq->len)) ||
       !_write(file, "\n", 1))
    {
      return 0;
    }

    len += seq->seq->len + 1;
  }
  else
  {
    unsigned long tmp_len = _write_wrapped(seq->seq, file);

    if(tmp_len == 0)
    {
      return 0;
    }

    len += tmp_len;
  }

  return len;
}

unsigned long _write_fastq_entry(SequenceFile *file, Sequence *seq)
{
  unsigned long len = 0;

  if(!_write(file, "@", 1) ||
     (seq->name->len != 0 && !_write(file, seq->name->buff, seq->name->len)) ||
     !_write(file, "\n", 1))
  {
    return 0;
  }

  len += 1 + seq->name->len + 1;

  // Check qual
  if(!seq->valid_qual)
  {
    unsigned long i;
    for(i = 0; i < seq->seq->len; i++)
    {
      string_buff_set_char(seq->qual, i, '?');
    }
  }

  if(file->line_wrap == 0)
  {
    if((seq->seq->len != 0 && !_write(file, seq->seq->buff, seq->seq->len)) ||
       !_write(file, "\n+\n", 3) ||
       (seq->qual->len != 0 && !_write(file, seq->qual->buff, seq->qual->len)) ||
       !_write(file, "\n", 1))
    {
      return 0;
    }

    len += seq->seq->len + 3 + seq->qual->len + 1;
  }
  else
  {
    unsigned long tmp_len;

    if(!(tmp_len = _write_wrapped(seq->seq, file)))
    {
      return 0;
    }

    len += tmp_len;

    if(!_write(file, "\n+\n", 3))
    {
      return 0;
    }

    len += 3;

    if(!(tmp_len = _write_wrapped(seq->qual, file)))
    {
      return 0;
    }

    len += tmp_len;
  }

  return len;
}

unsigned long _write_plain_entry(SequenceFile *file, Sequence *seq)
{
  if((seq->seq->len != 0 && !_write(file, seq->seq->buff, seq->seq->len)) ||
     !_write(file, "\n", 1))
  {
    return 0;
  }

  return seq->seq->len + 1;
}

unsigned long seq_file_write(SequenceFile *file, Sequence *seq)
{
  unsigned long chars_written;

  switch(file->file_type)
  {
    case SEQ_SAM:
    case SEQ_BAM:
      return 0;
    case SEQ_FASTA:
      chars_written = _write_fasta_entry(file, seq);
      break;
    case SEQ_FASTQ:
      chars_written = _write_fastq_entry(file, seq);
      break;
    case SEQ_PLAIN:
      chars_written = _write_plain_entry(file, seq);
      break;
    default:
      fprintf(stderr, "seq_reader.c warning: "
                      "tried to write to unknown filetype [path: %s]\n",
                      file->path);
      return 0;
  }

  if(chars_written > 0)
  {
    file->entry_number++;
  }

  return chars_written;
}

//
// Sequence Kmer Reader
//
SequenceKmerReader* kmer_reader_open(const char* path, const int kmer_size)
{
  SequenceFile* file = seq_file_open(path);

  if(file == NULL)
  {
    return NULL;
  }

  SequenceKmerReader* reader
    = (SequenceKmerReader*) malloc(sizeof(SequenceKmerReader));

  reader->file = file;
  reader->kmer_size = kmer_size;
  reader->offset = 0;

  if(file->file_type == SEQ_SAM ||
     file->file_type == SEQ_BAM ||
     file->file_type == SEQ_FASTQ)
  {
    reader->whole_entry = seq_create();
  }
  else
  {
    reader->whole_entry = NULL;
  }

  return reader;
}

void kmer_reader_close(SequenceKmerReader* reader)
{
  if(reader->whole_entry != NULL)
  {
    seq_free(reader->whole_entry);
  }

  free(reader);
}

long seq_file_entry_number(SequenceFile* file)
{
  return file->entry_number;
}

long kmer_reader_entry_number(SequenceKmerReader* reader)
{
  return reader->file->entry_number;
}

void _shift_insert_char(STRING_BUFFER* str, char c)
{
  memmove(str->buff, str->buff+1, str->len-1);
  str->buff[str->len-1] = c;
}

// Read an entry from a FASTA file
// Returns 1 if a header is read, 0 otherwise
char _read_fasta_entry_kmer(SequenceKmerReader* reader, Sequence* sequence)
{
  size_t kmer_size = reader->kmer_size;
  int c;

  if(reader->file->read_line_start == 1)
  {
    c = '>';
    reader->file->read_line_start = 0;
  }
  else
  {
    while((c = gzgetc(reader->file->gz_file)) != -1 && (c == '\n' || c == '\r'));

    if(c == -1)
    {
      seq_reset(sequence);
      return 0;
    }
  }

  if(string_buff_strlen(sequence->seq) == kmer_size && c != '>')
  {
    // Append
    _shift_insert_char(sequence->seq, (char)c);
    reader->offset++;
  }
  else
  {
    // Read new
    seq_reset(sequence);

    if(c != '>')
    {
      fprintf(stderr, "seq_reader.c - read_fasta_entry_kmer(): expected "
                      "line to start '>'\n");
      seq_reset(sequence);
      return 0;
    }

    while(string_buff_strlen(sequence->seq) < kmer_size && c != -1)
    {
      // reset
      seq_reset(sequence);
      reader->offset = 0;
      reader->file->entry_number++;

      // Read name
      string_buff_gzreadline(sequence->name, reader->file->gz_file);
      string_buff_chomp(sequence->name);

      // Read sequence over multiple lines
      while(string_buff_strlen(sequence->seq) < kmer_size)
      {
        while((c = gzgetc(reader->file->gz_file)) != -1 &&
              (c == '\n' || c == '\r'));

        if(c == -1 || c == '>')
        {
          break;
        }

        string_buff_append_char(sequence->seq, c);
        string_buff_gzread(sequence->seq, reader->file->gz_file,
                           kmer_size - string_buff_strlen(sequence->seq));

        string_buff_chomp(sequence->seq);
      }

      //printf(">%s:%c:'%s'\n", sequence->name->buff, c, sequence->seq->buff);
    }

    if(c == '>')
    {
      // No sequence with name, but still success
      reader->file->read_line_start = 1;
    }

    if(string_buff_strlen(sequence->seq) < kmer_size)
    {
      seq_reset(sequence);
      return 0;
    }
  }

  sequence->offset = reader->offset;

  return 1;
}

// Read an entry from a SAM/BAM file
// Returns 1 if a line is read, 0 otherwise
char _read_plain_entry_kmer(SequenceKmerReader* reader, Sequence* sequence)
{
  // Plain is one per sequence
  // Read a single char
  int c;
  char new_line = 0;

  if(reader->file->read_line_start)
  {
    // First entry only
    c = reader->file->read_line_start;
    reader->file->read_line_start = 0;
    new_line = 1;
  }
  else if((c = gzgetc(reader->file->gz_file)) == -1)
  {
    seq_reset(sequence);
    return 0;
  }

  if(!new_line && string_buff_strlen(sequence->seq) == reader->kmer_size &&
     c != '\n' && c != '\r')
  {
    // Append
    _shift_insert_char(sequence->seq, (char)c);
    reader->offset++;
  }
  else
  {
    // Read new
    seq_reset(sequence);
    reader->offset = 0;
    reader->file->entry_number++;

    if(c != '\n' && c != '\r')
    {
      // Append char
      string_buff_append_char(sequence->seq, c);
    }

    t_buf_pos remaining;

    while(1)
    {
      remaining = reader->kmer_size - string_buff_strlen(sequence->seq);

      if(string_buff_gzread(sequence->seq,
                            reader->file->gz_file,
                            remaining) == 0)
      {
        // Hit end of file
        seq_reset(sequence);
        return 0;
      }

      string_buff_chomp(sequence->seq);

      if(string_buff_strlen(sequence->seq) < reader->kmer_size)
      {
        seq_reset(sequence);
        reader->file->entry_number++;
      }
      else
      {
        break;
      }
    }
  }

  sequence->offset = reader->offset;
  sequence->valid_qual = 0;

  return 1;
}

// Read SAM/BAM or FASTQ
char _read_whole_entry_kmer(SequenceKmerReader *reader, Sequence *sequence)
{
  // Need to read in a whole sequence for these ones
  char read_new_entry = 0;

  if(string_buff_strlen(reader->whole_entry->seq) > 0)
  {
    // If we've read in a sequence, advance the offset
    reader->offset++;
  }

  while(reader->offset + reader->kmer_size >
        string_buff_strlen(reader->whole_entry->seq))
  {
    // Ran out of sequence on current entry
    seq_reset(reader->whole_entry);
    reader->offset = 0;
    reader->file->entry_number++;

    if(reader->file->file_type == SEQ_FASTQ)
    {
      read_new_entry = _read_fastq_entry(reader->file, reader->whole_entry);
    }
    else
    {
      // SAM or BAM
      read_new_entry = _read_bam_entry(reader->file, reader->whole_entry);
    }

    if(!read_new_entry)
    {
      seq_reset(sequence);
      return 0;
    }
  }

  if(read_new_entry)
  {
    seq_reset(sequence);

    // copy name
    string_buff_copy(sequence->name, 0, reader->whole_entry->name, 0,
                     string_buff_strlen(reader->whole_entry->name));
    // copy sequence
    string_buff_copy(sequence->seq, 0, reader->whole_entry->seq, 0,
                     reader->kmer_size);
    // copy quality
    if(string_buff_strlen(reader->whole_entry->qual) > 0)
    {
      string_buff_copy(sequence->qual, 0, reader->whole_entry->qual, 0,
                       reader->kmer_size);
    }
    return 1;
  }
  else
  {
    // append
    t_buf_pos pos = reader->offset + reader->kmer_size - 1;
    char next;

    next = string_buff_get_char(reader->whole_entry->seq, pos);
    _shift_insert_char(sequence->seq, next);

    next = string_buff_get_char(reader->whole_entry->qual, pos);
    _shift_insert_char(sequence->qual, next);

    sequence->offset = reader->offset;
    return 1;
  }
}

// Returns 0 if at end of file; 1 otherwise
char kmer_reader_read_kmer(SequenceKmerReader* reader, Sequence* sequence)
{
  SequenceFileType file_type = reader->file->file_type;

  if(file_type == SEQ_SAM || file_type == SEQ_BAM || file_type == SEQ_FASTQ)
  {
    return _read_whole_entry_kmer(reader, sequence);
  }
  else if(file_type == SEQ_FASTA)
  {
    return _read_fasta_entry_kmer(reader, sequence);
  }
  else if(file_type == SEQ_PLAIN)
  {
    return _read_plain_entry_kmer(reader, sequence);
  }
  else
  {
    fprintf(stderr, "seq_reader.c warning: "
                      "tried to read kmer from unknown filetype\n");
    return 0;
  }
}

int _read_whole_entry_char(SequenceKmerReader* reader)
{
  if(reader->offset < string_buff_strlen(reader->whole_entry->seq))
  {
    return string_buff_get_char(reader->whole_entry->seq, reader->offset);
  }
  else
  {
    return -1;
  }
}

int _read_fasta_entry_char(SequenceKmerReader* reader)
{
  if(reader->file->read_line_start)
  {
    return -1;
  }

  char c;
  
  while((c = gzgetc(reader->file->gz_file)) != -1 && (c == '\n' || c == '\r'));

  if(c == -1)
  {
    return -1;
  }
  else if(c == '>')
  {
    reader->file->read_line_start = 1;
    return -1;
  }

  return c;
}

int _read_plain_entry_char(SequenceKmerReader* reader)
{
  if(reader->file->read_line_start)
  {
    return -1;
  }

  char c = gzgetc(reader->file->gz_file);

  if(c == '\n' || c == '\r')
  {
    reader->file->read_line_start = 1;
    return -1;
  }

  // will return -1 if we hit the end
  return c;
}

// Return next char or -1 if end of read entry
int kmer_reader_read_char(SequenceKmerReader* reader)
{
  SequenceFileType file_type = reader->file->file_type;

  if(reader->offset == 0)
  {
    reader->offset = reader->kmer_size;
  }
  else
  {
    reader->offset++;
  }

  if(file_type == SEQ_SAM || file_type == SEQ_BAM || file_type == SEQ_FASTQ)
  {
    return _read_whole_entry_char(reader);
  }
  else if(file_type == SEQ_FASTA)
  {
    return _read_fasta_entry_char(reader);
  }
  else if(file_type == SEQ_PLAIN)
  {
    return _read_plain_entry_char(reader);
  }
  else
  {
    fprintf(stderr, "seq_reader.c warning: "
                    "tried to read kmer from unknown filetype [path: %s]\n",
                    reader->file->path);
    return 0;
  }
}
