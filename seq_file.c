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

#include "seq_file.h"
#include "sam.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))

#define is_base_char(x) ((x) == 'a' || (x) == 'A' || \
                         (x) == 'c' || (x) == 'C' || \
                         (x) == 'g' || (x) == 'G' || \
                         (x) == 't' || (x) == 'T' || \
                         (x) == 'n' || (x) == 'N')

// Write output MACROs
// wrapper for fputs/gzputs
#define seq_puts(f,str) (unsigned long) \
((f)->plain_file != NULL ? fputs((str), (f)->plain_file) \
                         : gzputs((f)->gz_file, (str)))

// wrapper for fwrite/gzwrite
#define seq_write(f,str,len) \
((f)->plain_file != NULL \
  ? (unsigned long)fwrite((str), sizeof(char), (len), (f)->plain_file) \
  : (unsigned long)gzwrite((f)->gz_file, (str), (len)))

const char* seq_file_types[6]
  = {"Unknown", "FASTA", "FASTQ", "Plain", "SAM", "BAM"};

const char* seq_file_types_zipped[6]
  = {"Unknown(zipped)", "FASTA(zipped)", "FASTQ(zipped)", "Plain(zipped)",
     "SAM", "BAM"};

// Array for complementing bases read from BAM/SAM files
int8_t seq_comp_table[16] = {0,8,4,12,2,10,9,14,1,6,5,13,3,11,7,15};

// Printed nothing, then name, then some sequence, then some qualities
// ... then ready to print a name again
typedef enum WriteState
  {WS_READ_ONLY, WS_BEGIN, WS_NAME, WS_SEQ, WS_QUAL} WriteState;

struct SeqFile
{
  const char *path;

  gzFile *gz_file; // for reading FASTA/FASTQ/plain
  samfile_t *sam_file; // For reading SAM/BAM

  // For reading sam/bams
  bam1_t *bam;

  enum SeqFileType file_type;

  // have we seen a '>' at the start of a line in a fasta file?
  char read_line_start;

  // name, index and bases-read/offset of current entry
  STRING_BUFFER *entry_name;
  unsigned long entry_index;
  
  unsigned long entry_offset, entry_offset_qual;

  // Whether an entry has been read in
  char entry_read, entry_read_qual;

  // Buffer for reading in bases in FASTQ files
  STRING_BUFFER *bases_buff;

  // total bases read/written - initially 0
  unsigned long total_bases_passed;

  /* Writing preferences */

  // Output plain file for writing if not gzipping output
  // (newer zlib allows you to do this with gzFile, but mac version is outdated)
  FILE *plain_file;

  // 0 if no wrap, otherwise max bases per line
  unsigned long line_wrap, curr_line_length;

  // State of writing
  WriteState write_state;
};

//
// Sequence File reader
//

void seq_guess_filetype_from_path(const char *path, SeqFileType *file_type,
                                  char *zipped)
{
  size_t path_len = strlen(path);

  #define NUM_FA 2
  #define NUM_GZ 5

  const char* fa[NUM_FA] = {".fa",".fasta"};
  const char* fq[NUM_FA] = {".fq",".fastq"};
  const char* fagz[NUM_GZ] = {".faz",".fagz",".fa.gz",".fa.gzip",".fasta.gzip"};
  const char* fqgz[NUM_GZ] = {".fqz",".fqgz",".fq.gz",".fq.gzip",".fastq.gzip"};

  if(strcasecmp(".sam", path+path_len - 4) == 0)
  {
    *file_type = SEQ_SAM;
    return;
  }
  else if(strcasecmp(".bam", path+path_len - 4) == 0)
  {
    *file_type = SEQ_BAM;
    return;
  }
  
  int i;

  // FASTA
  for(i = 0; i < NUM_FA; i++)
  {
    if(strcasecmp(fa[i], path+path_len-strlen(fa[i])) == 0)
    {
      *zipped = 0;
      *file_type = SEQ_FASTA;
      return;
    }
  }

  // FASTA ZIPPED
  for(i = 0; i < NUM_GZ; i++)
  {
    if(strcasecmp(fagz[i], path+path_len-strlen(fagz[i])) == 0)
    {
      *zipped = 1;
      *file_type = SEQ_FASTA;
      return;
    }
  }
  
  // FASTA
  for(i = 0; i < NUM_FA; i++)
  {
    if(strcasecmp(fq[i], path+path_len-strlen(fq[i])) == 0)
    {
      *zipped = 0;
      *file_type = SEQ_FASTQ;
      return;
    }
  }

  // FASTA ZIPPED
  for(i = 0; i < NUM_GZ; i++)
  {
    if(strcasecmp(fqgz[i], path+path_len-strlen(fqgz[i])) == 0)
    {
      *zipped = 1;
      *file_type = SEQ_FASTQ;
      return;
    }
  }

  if(strcasecmp(".txt", path+path_len-4) == 0)
  {
    *zipped = 0;
    *file_type = SEQ_PLAIN;
    return;
  }

  if(strcasecmp(".txt.gz", path+path_len-7) == 0)
  {
    *zipped = 1;
    *file_type = SEQ_PLAIN;
    return;
  }

  *file_type = SEQ_UNKNOWN;
}

void _set_seq_filetype(SeqFile *file)
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

  if(file->gz_file == NULL)
  {
    return;
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
    file->bases_buff = string_buff_new();
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

SeqFile* _create_default_seq_file(const char* file_path)
{
  SeqFile* sf = (SeqFile*) malloc(sizeof(SeqFile));

  sf->path = file_path;

  sf->gz_file = NULL;
  sf->sam_file = NULL;

  sf->bam = NULL;

  sf->file_type = SEQ_UNKNOWN;
  sf->read_line_start = 0;

  sf->entry_name = string_buff_new();
  sf->entry_index = 0;

  sf->entry_offset = 0;
  sf->entry_offset_qual = 0;

  sf->entry_read = 0;
  sf->entry_read_qual = 0;

  sf->bases_buff = NULL;

  sf->total_bases_passed = 0;

  // For writing
  sf->plain_file = NULL;
  sf->line_wrap = 0;
  sf->curr_line_length = 0;
  sf->write_state = WS_READ_ONLY;

  return sf;
}

SeqFile* seq_file_open(const char* file_path)
{
  SeqFile* file = _create_default_seq_file(file_path);
  _set_seq_filetype(file);

  if(file->file_type == SEQ_UNKNOWN)
  {
    free(file);
    return NULL;
  }
  else
  {
    return file;
  }
}

SeqFile* seq_file_open_filetype(const char* file_path,
                                const SeqFileType file_type)
{
  SeqFile* file = _create_default_seq_file(file_path);
  file->file_type = file_type;

  if(file_type == SEQ_FASTQ)
  {
    file->bases_buff = string_buff_new();
  }

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
      fprintf(stderr, "seq_reader.c warning: invalid SeqFileType in "
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
SeqFile* seq_file_open_write(const char* file_path, const SeqFileType file_type,
                             const char gzip, unsigned long line_wrap)
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
    line_wrap = 0;
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

  SeqFile* sf = _create_default_seq_file(file_path);

  sf->gz_file = gz_file;
  sf->plain_file = plain_file;
  sf->file_type = file_type;
  sf->line_wrap = line_wrap;

  // Set write state (default is read-only)
  sf->write_state = WS_BEGIN;

  return sf;
}


// Close
// Adds an extra byte (newline) to the end of output (writing) files
size_t seq_file_close(SeqFile* sf)
{
  size_t num_bytes_printed = 0;

  // Add a new line to the end of output files
  if(sf->write_state != WS_READ_ONLY)
  {
    num_bytes_printed = seq_puts(sf, "\n");
  }

  if(sf->gz_file != NULL)
  {
    gzclose(sf->gz_file);
  }

  if(sf->file_type == SEQ_SAM || sf->file_type == SEQ_BAM)
  {
    bam_destroy1(sf->bam);
    samclose(sf->sam_file);
  }

  if(sf->bases_buff != NULL)
  {
    free(sf->bases_buff);
  }

  if(sf->plain_file != NULL)
  {
    fclose(sf->plain_file);
  }

  string_buff_free(sf->entry_name);

  free(sf);

  return num_bytes_printed;
}

SeqFileType seq_file_get_type(const SeqFile* file)
{
  return file->file_type;
}

const char* seq_file_get_type_str(const SeqFile* file)
{
  return seq_file_types[file->file_type];
}

const char* seq_file_type_str(const SeqFileType file_type,
                              const char zipped)
{
  return zipped ? seq_file_types_zipped[file_type] : seq_file_types[file_type];
}

// Get a pointer to the file path
const char* seq_get_path(const SeqFile* file)
{
  return file->path;
}

unsigned long seq_num_bases_passed(const SeqFile *sf)
{
  return sf->total_bases_passed;
}

char seq_has_quality_scores(const SeqFile *sf)
{
  switch (sf->file_type)
  {
    case SEQ_SAM:
    case SEQ_BAM:
    case SEQ_FASTQ:
      return 1;
    case SEQ_FASTA:
    case SEQ_PLAIN:
      return 0;
    default:
      fprintf(stderr, "seq_file.c: unknown file type [path: %s]\n", sf->path);
      return 0;
  }
}

// Returns 1 if open for writing, 0 otherwise
char seq_is_open_for_write(const SeqFile *sf)
{
  return (sf->write_state != WS_READ_ONLY);
}

char _seq_next_read_fasta(SeqFile *sf)
{
  if(sf->read_line_start)
  {
    // Read name
    string_buff_gzreadline(sf->entry_name, sf->gz_file);
    string_buff_chomp(sf->entry_name);

    sf->read_line_start = 0;
    return 1;
  }
  else
  {
    int c;

    // Look for line starting with >
    do
    {
      // Read until the end of the line
      while((c = gzgetc(sf->gz_file)) != -1 && c != '\n' && c != '\r');

      if(c == -1)
        return 0;

      // Read through end of line chars
      while((c = gzgetc(sf->gz_file)) != -1 && (c == '\n' || c == '\r'));

      if(c == -1)
        return 0;
    }
    while(c != '>');

    // Read name
    string_buff_gzreadline(sf->entry_name, sf->gz_file);
    string_buff_chomp(sf->entry_name);

    return 1;
  }
}

void _seq_read_fastq_sequence(SeqFile *sf)
{
  string_buff_reset(sf->bases_buff);

  int c;

  while((c = gzgetc(sf->gz_file)) != -1 && c != '+')
  {
    if(c != '\r' && c != '\n')
    {
      string_buff_append_char(sf->bases_buff, c);
      string_buff_gzreadline(sf->bases_buff, sf->gz_file);
      string_buff_chomp(sf->bases_buff);
    }
  }

  if(c == -1)
  {
    fprintf(stderr, "seq_file.c: missing + in FASTQ [file: %s]\n", sf->path);
  }
}

char _seq_next_read_fastq(SeqFile *sf)
{
  if(sf->read_line_start)
  {
    // Read name
    string_buff_gzreadline(sf->entry_name, sf->gz_file);
    string_buff_chomp(sf->entry_name);

    // Read whole sequence
    _seq_read_fastq_sequence(sf);

    sf->read_line_start = 0;
    return 1;
  }
  else
  {
    int c;

    // Skip over remaining quality values
    while(sf->entry_offset_qual < string_buff_strlen(sf->bases_buff))
    {
      if((c = gzgetc(sf->gz_file)) == -1)
        return 0;

      if(c != '\r' && c != '\n')
        sf->entry_offset_qual++;
    }

    // Skip newlines
    while((c = gzgetc(sf->gz_file)) != -1 && (c == '\n' || c == '\r'));

    if(c == -1)
      return 0;

    if(c != '@')
    {
      fprintf(stderr, "seq_file.c: FASTQ header does not begin with '@' [%c]\n",
              c);
      return 0;
    }

    // Read name
    string_buff_gzreadline(sf->entry_name, sf->gz_file);
    string_buff_chomp(sf->entry_name);

    // Read whole sequence
    _seq_read_fastq_sequence(sf);

    return 1;
  }
}

char _seq_next_read_bam(SeqFile *sf)
{
  if(samread(sf->sam_file, sf->bam) < 0)
    return 0;

  // Get name
  string_buff_append_str(sf->entry_name, bam1_qname(sf->bam));

  return 1;
}

char _seq_next_read_plain(SeqFile *sf)
{
  if(sf->read_line_start)
  {
    sf->read_line_start = 0;
    return 1;
  }
  else
  {
    int c;

    while((c = gzgetc(sf->gz_file)) != -1 && c != '\r' && c != '\n');

    if(c == -1)
      return 0;

    return 1;
  }
}

// Returns 1 on success 0 if no more to read
char seq_next_read(SeqFile *sf)
{
  string_buff_reset(sf->entry_name);

  char success;

  switch (sf->file_type)
  {
    case SEQ_FASTA:
      success = _seq_next_read_fasta(sf);
      break;
    case SEQ_FASTQ:
      success = _seq_next_read_fastq(sf);
      break;
    case SEQ_SAM:
    case SEQ_BAM:
      success = _seq_next_read_bam(sf);
      break;
    case SEQ_PLAIN:
      success = _seq_next_read_plain(sf);
      break;
    default:
      fprintf(stderr, "seq_file.c: Cannot read from unknown file type "
                      "[file: %s]\n", sf->path);
      return 0;
  }

  sf->entry_offset = 0;
  sf->entry_offset_qual = 0;

  if(success)
  {
    sf->entry_index++;
    sf->entry_read = 1;
    sf->entry_read_qual = seq_has_quality_scores(sf);
  }
  else
  {
    sf->entry_read = 0;
    sf->entry_read_qual = 0;
  }

  return success;
}


// Get the name of the next read
const char* seq_get_read_name(SeqFile *sf)
{
  return sf->entry_name->buff;
}

// Get this read index -- starts from 0
unsigned long seq_get_read_index(SeqFile *sf)
{
  return sf->entry_index;
}

unsigned long seq_get_base_offset(SeqFile *sf)
{
  return sf->entry_offset;
}

unsigned long seq_get_qual_offset(SeqFile *sf)
{
  return sf->entry_offset_qual;
}

unsigned long seq_get_length(SeqFile *sf)
{
  if(sf->entry_read)
  {
    fprintf(stderr, "seq_file.c: haven't finished reading sequence - "
                    "seq_get_length() cannot return a length\n");
    return 0;
  }

  return sf->entry_offset;
}


/*
 Read a single base from a read
*/

char _seq_read_base_fasta(SeqFile *sf, char *c)
{
  int next;
  
  while((next = gzgetc(sf->gz_file)) != -1 && (next == '\n' || next == '\r'));

  if(next == -1)
  {
    return 0;
  }
  else if(next == '>')
  {
    sf->read_line_start = 1;
    return 0;
  }
  else
  {
    *c = next;
    return 1;
  }
}

char _seq_read_base_fastq(SeqFile *sf, char *c)
{
  if(sf->entry_offset < string_buff_strlen(sf->bases_buff))
  {
    *c = string_buff_get_char(sf->bases_buff, sf->entry_offset);
    return 1;
  }
  else
  {
    return 0;
  }
}

char _seq_read_base_bam(SeqFile *sf, char *c)
{
  // Get reverse
  char is_reversed = sf->bam->core.flag & 16;
  int query_len = sf->bam->core.l_qseq;

  if(sf->entry_offset >= (unsigned long)query_len)
    return 0;

  uint8_t *seq = bam1_seq(sf->bam);
  
  if(is_reversed)
  {
    int index = query_len - sf->entry_offset - 1;
    int8_t b = bam1_seqi(seq, index);
    *c = bam_nt16_rev_table[seq_comp_table[b]];
  }
  else
  {
    int8_t b = bam1_seqi(seq, sf->entry_offset);
    *c = bam_nt16_rev_table[b];
  }

  return 1;
}

char _seq_read_base_plain(SeqFile *sf, char *c)
{
  int next = gzgetc(sf->gz_file);

  if(next == -1)
  {
    return 0;
  }
  else if(next == '\r' || next == '\n')
  {
    sf->read_line_start = 1;
    return 0;
  }
  else
  {
    *c = next;
    return 1;
  }
}

// Read a single base from the current read
// Returns 1 on success, 0 if no more quality scores or run out of bases
char seq_read_base(SeqFile *sf, char *c)
{
  // Check if we have read anything in
  if(!sf->entry_read)
    return 0;

  char success;

  switch (sf->file_type)
  {
    case SEQ_FASTA:
      success = _seq_read_base_fasta(sf, c);
      break;
    case SEQ_FASTQ:
      success = _seq_read_base_fastq(sf, c);
      break;
    case SEQ_SAM:
    case SEQ_BAM:
      success = _seq_read_base_bam(sf, c);
      break;
    case SEQ_PLAIN:
      success = _seq_read_base_plain(sf, c);
      break;
    default:
      fprintf(stderr, "seq_file.c: Cannot read from unknown file type "
                      "[file: %s]\n", sf->path);
      return 0;
  }

  if(success)
  {
    sf->entry_offset++;
    sf->total_bases_passed++;
  }

  sf->entry_read = success;

  return success;
}


/*
 Read a single quality score from a read
*/

char _seq_read_qual_fastq(SeqFile *sf, char *c)
{
  if(sf->entry_offset_qual >= string_buff_strlen(sf->bases_buff))
    return 0;

  int next;

  while((next = gzgetc(sf->gz_file)) != -1 && (next == '\n' || next == '\r'));

  if(next == -1)
  {
    fprintf(stderr, "seq_file.c: fastq file ended without finishing quality "
                    "scores [file: %s]\n", sf->path);
    return 0;
  }

  *c = next;
  return 1;
}

char _seq_read_qual_bam(SeqFile *sf, char *c)
{
  char is_reversed = sf->bam->core.flag & 16;
  int query_len = sf->bam->core.l_qseq;

  if(sf->entry_offset_qual >= (unsigned long)query_len)
    return 0;

  uint8_t *seq = bam1_qual(sf->bam);
  int index;

  if(is_reversed)
    index = query_len - sf->entry_offset_qual - 1;
  else
    index = sf->entry_offset_qual;

  *c = 33 + seq[index];

  return 1;
}

// Read a single quality score from the current read
// Returns 1 on success, 0 if no more quality scores or run out of bases
char seq_read_qual(SeqFile *sf, char *c)
{
  // Check if we have read anything in
  if(!sf->entry_read_qual)
    return 0;

  char success;

  switch (sf->file_type)
  {
    case SEQ_FASTQ:
      success = _seq_read_qual_fastq(sf, c);
      break;
    case SEQ_SAM:
    case SEQ_BAM:
      success = _seq_read_qual_bam(sf, c);
      break;
    case SEQ_FASTA:
    case SEQ_PLAIN:
      fprintf(stderr, "seq_file.c: Cannot read from file type without quality "
                      "scores [file: %s; type: %s]\n", sf->path,
                      seq_file_type_str(sf->file_type, 0));
      return 0;
    default:
      fprintf(stderr, "seq_file.c: Cannot read from unknown file type "
                      "[file: %s]\n", sf->path);
      return 0;
  }

  if(success)
    sf->entry_offset_qual++;

  sf->entry_read_qual = success;

  return success;
}

/*
 Read k bases / quality scores from a read
*/

// str must be at least k+1 bytes long
// returns 1 on success, 0 otherwise
char seq_read_k_bases(SeqFile *sf, char* str, int k)
{
  // Check if we have read anything in
  if(!sf->entry_read)
    return 0;

  // read will have been exhausted once this method returns
  sf->entry_read = 0;

  int i;

  for(i = 0; i < k; i++)
  {
    if(!seq_read_base(sf, str+i))
    {
      str[0] = '\0';
      return 0;
    }
  }

  str[k] = '\0';

  return 1;
}

char seq_read_k_quals(SeqFile *sf, char* str, int k)
{
  // Check if we have read anything in
  if(!sf->entry_read_qual)
    return 0;

  // read will have been exhausted once this method returns
  sf->entry_read_qual = 0;

  int i;

  for(i = 0; i < k; i++)
  {
    if(!seq_read_qual(sf, str+i))
    {
      str[0] = '\0';
      return 0;
    }
  }

  str[k] = '\0';
  return 1;
}


/*
 Read the rest of a read
*/

char _seq_read_all_bases_bam(SeqFile *sf, STRING_BUFFER *sbuf)
{
  int qlen = sf->bam->core.l_qseq;

  if(sf->entry_offset >= (unsigned long)qlen)
    return 0;

  // Get reverse
  char is_reversed = sf->bam->core.flag & 16;

  string_buff_ensure_capacity(sbuf, qlen - sf->entry_offset);

  uint8_t *seq = bam1_seq(sf->bam);

  // read in and reverse complement (if needed)
  int i;
  for(i = sf->entry_offset; i < qlen; i++)
  {
    int index = (is_reversed ? i : qlen - i - 1);
    int8_t b = bam1_seqi(seq, index);
    char c = bam_nt16_rev_table[is_reversed ? seq_comp_table[b] : b];
    string_buff_append_char(sbuf, c);
  }

  return 1;
}

char _seq_read_all_bases_fasta(SeqFile *sf, STRING_BUFFER *sbuf)
{
  int c;

  while((c = gzgetc(sf->gz_file)) != -1 && c != '>')
  {
    if(c != '\r' && c != '\n')
    {
      string_buff_append_char(sbuf, c);
      string_buff_gzreadline(sbuf, sf->gz_file);
      string_buff_chomp(sbuf);
    }
  }

  if(c == '>')
    sf->read_line_start = 1;

  return 1;
}

char _seq_read_all_bases_fastq(SeqFile *sf, STRING_BUFFER *sbuf)
{
  // Copy from buffer
  t_buf_pos len = sf->bases_buff->len - sf->entry_offset;
  string_buff_copy(sbuf, 0, sf->bases_buff, sf->entry_offset, len);

  return 1;
}

char _seq_read_all_bases_plain(SeqFile *sf, STRING_BUFFER *sbuf)
{
  t_buf_pos len = string_buff_gzreadline(sbuf, sf->gz_file);
  string_buff_chomp(sbuf);

  return (len > 0);
}

// returns 1 on success, 0 otherwise
char seq_read_all_bases(SeqFile *sf, STRING_BUFFER *sbuf)
{
  // Check if we have read anything in
  if(!sf->entry_read)
    return 0;

  string_buff_reset(sbuf);

  char success;

  switch(sf->file_type)
  {
    case SEQ_SAM:
    case SEQ_BAM:
      success = _seq_read_all_bases_bam(sf, sbuf);
      break;
    case SEQ_FASTA:
      success = _seq_read_all_bases_fasta(sf, sbuf);
      break;
    case SEQ_FASTQ:
      success = _seq_read_all_bases_fastq(sf, sbuf);
      break;
    case SEQ_PLAIN:
      success = _seq_read_all_bases_plain(sf, sbuf);
      break;
    default:
      fprintf(stderr, "seq_file.c: tried to read from unknown filetype "
                      "[path: %s]\n", sf->path);
      return 0;
  }

  sf->entry_offset += string_buff_strlen(sbuf);
  sf->total_bases_passed += string_buff_strlen(sbuf);

  // read has been exhausted
  sf->entry_read = 0;

  return success;
}


/*
 Read the rest of a read's quality scores
*/

char _seq_read_all_quals_bam(SeqFile *sf, STRING_BUFFER *sbuf)
{
  int qlen = sf->bam->core.l_qseq;

  if(sf->entry_offset_qual >= (unsigned long)qlen)
    return 0;

  // Get reverse
  char is_reversed = sf->bam->core.flag & 16;

  string_buff_ensure_capacity(sbuf, qlen - sf->entry_offset);

  uint8_t *seq = bam1_qual(sf->bam);

  // read in and reverse complement (if needed)
  int i;
  for(i = sf->entry_offset; i < qlen; i++)
  {
    char c = 33 + seq[is_reversed ? i : qlen - i - 1];
    string_buff_append_char(sbuf, c);
  }

  return 1;
}

char _seq_read_all_quals_fastq(SeqFile *sf, STRING_BUFFER *sbuf)
{
  if(sf->entry_offset_qual >= string_buff_strlen(sf->bases_buff))
    return 0;

  // Expect the same number of quality scores as bases
  t_buf_pos expected_len = string_buff_strlen(sf->bases_buff) -
                           sf->entry_offset_qual;

  int next = -1;
  t_buf_pos i;

  for(i = 0; i < expected_len && (next = gzgetc(sf->gz_file)) != -1; i++)
  {
    if(next != '\r' && next != '\n')
    {
      string_buff_append_char(sbuf, next);
    }
  }

  if(next == -1)
  {
    fprintf(stderr, "seq_file.c: fastq file ended without finishing quality "
                    "scores (FASTQ) [file: %s]\n", sf->path);
  }

  return 1;
}

// returns 1 on success, 0 otherwise
char seq_read_all_quals(SeqFile *sf, STRING_BUFFER *sbuf)
{
  // Check if we have read anything in
  if(!sf->entry_read_qual)
    return 0;

  string_buff_reset(sbuf);

  char success;

  switch(sf->file_type)
  {
    case SEQ_SAM:
    case SEQ_BAM:
      success = _seq_read_all_quals_bam(sf, sbuf);
      break;
    case SEQ_FASTQ:
      success = _seq_read_all_quals_fastq(sf, sbuf);
      break;
  default:
      fprintf(stderr, "seq_file.c: tried to read from unknown filetype "
                      "[path: %s]\n", sf->path);
      return 0;
  }

  // Exhausted read quality scores
  sf->entry_offset_qual += string_buff_strlen(sbuf);
  sf->entry_read_qual = 0;

  return success;
}


/*
 Write to a file. 
 Each function returns the number of bytes written or 0 on failure
*/

unsigned long _seq_file_write_name_fasta(SeqFile *sf, const char *name)
{
  size_t num_bytes_printed = 0;

  if(sf->write_state == WS_BEGIN)
  {
    num_bytes_printed += seq_puts(sf, ">");
  }
  else if(sf->write_state == WS_SEQ)
  {
    num_bytes_printed += seq_puts(sf, "\n>");
  }

  num_bytes_printed += seq_puts(sf, name);

  return num_bytes_printed;
}

unsigned long _seq_file_write_name_fastq(SeqFile *sf, const char *name)
{
  if(sf->write_state == WS_SEQ)
  {
    fprintf(stderr, "seq_file.c: writing in the wrong order (name) [path: %s]\n",
            sf->path);
    exit(EXIT_FAILURE);
  }

  unsigned long num_bytes_printed = 0;

  if(sf->write_state == WS_BEGIN)
  {
    num_bytes_printed += seq_puts(sf, "@");
  }
  else if(sf->write_state == WS_QUAL)
  {
    num_bytes_printed += seq_puts(sf, "\n@");
  }

  num_bytes_printed += seq_puts(sf, name);

  return num_bytes_printed;
}

unsigned long seq_file_write_name(SeqFile *sf, const char *name)
{
  unsigned long num_bytes_printed = 0;

  switch(sf->file_type)
  {
    case SEQ_FASTA:
      num_bytes_printed = _seq_file_write_name_fasta(sf, name);
      break;
    case SEQ_FASTQ:
      num_bytes_printed = _seq_file_write_name_fastq(sf, name);
      break;
    default:
      fprintf(stderr, "seq_file.c: called seq_file_write_name() with invalid "
                      "file type (%s) [path: %s]\n", seq_file_get_type_str(sf),
                      sf->path);
      return 0;
  }

  sf->write_state = WS_NAME;
  return num_bytes_printed;
}



/*
 Write sequence
*/

#define _write(s,str,len) \
((sf)->line_wrap == 0 ? seq_puts((sf), (str)) : _write_wrapped((sf),(str),(len)))

inline size_t _write_wrapped(SeqFile *sf, const char *str, size_t str_len)
{
  size_t num_bytes_printed = 0;

  if(sf->curr_line_length == sf->line_wrap)
  {
    sf->curr_line_length = 0;
    num_bytes_printed += seq_puts(sf, "\n");
  }
  
  if(sf->curr_line_length + str_len <= sf->line_wrap)
  {
    // Doesn't go over a single line
    sf->curr_line_length += str_len;
    num_bytes_printed += seq_puts(sf, str);
    return num_bytes_printed;
  }

  size_t bytes_to_print = sf->line_wrap - sf->curr_line_length;

  num_bytes_printed += seq_write(sf, str, bytes_to_print);
  num_bytes_printed += seq_puts(sf, "\n");

  size_t offset;

  for(offset = bytes_to_print; offset < str_len; offset += sf->line_wrap)
  {
    bytes_to_print = MIN(str_len - offset, sf->line_wrap);
    num_bytes_printed += seq_write(sf, str + offset, bytes_to_print);

    if(bytes_to_print < sf->line_wrap)
      num_bytes_printed += seq_puts(sf, "\n");
  }

  sf->curr_line_length = bytes_to_print;

  return num_bytes_printed;
}

// Write FASTA sequence
size_t _seq_file_write_seq_fasta(SeqFile *sf, const char *seq, size_t str_len)
{
  if(sf->write_state == WS_BEGIN)
  {
    fprintf(stderr, "seq_file.c: writing in the wrong order (seq) [path: %s]\n",
            sf->path);
    exit(EXIT_FAILURE);
  }

  size_t num_bytes_printed = 0;

  if(sf->write_state == WS_NAME)
  {
    num_bytes_printed += seq_puts(sf, "\n");
  }

  num_bytes_printed += _write(sf, seq, str_len);

  return num_bytes_printed;
}

size_t _seq_file_write_seq_fastq(SeqFile *sf, const char *seq, size_t str_len)
{
  if(sf->write_state == WS_BEGIN || sf->write_state == WS_QUAL)
  {
    fprintf(stderr, "seq_file.c: writing in the wrong order (seq) [path: %s]\n",
            sf->path);
    exit(EXIT_FAILURE);
  }

  size_t num_bytes_printed = 0;

  if(sf->write_state == WS_NAME)
  {
    num_bytes_printed += seq_puts(sf, "\n");
  }

  num_bytes_printed += _write(sf, seq, str_len);

  return num_bytes_printed;
}

size_t _seq_file_write_seq_plain(SeqFile *sf, const char *seq)
{
  size_t num_bytes_printed = 0;

  if(sf->write_state != WS_BEGIN)
    num_bytes_printed += seq_puts(sf, "\n");

  num_bytes_printed += seq_puts(sf, seq);

  return num_bytes_printed;
}

size_t seq_file_write_seq(SeqFile *sf, const char *seq)
{
  size_t str_len = strlen(seq);

  if(str_len == 0)
    return 0;

  unsigned long num_bytes_printed = 0;

  switch(sf->file_type)
  {
    case SEQ_FASTA:
      num_bytes_printed = _seq_file_write_seq_fasta(sf, seq, str_len);
      break;
    case SEQ_FASTQ:
      num_bytes_printed = _seq_file_write_seq_fastq(sf, seq, str_len);
      break;
    case SEQ_PLAIN:
      num_bytes_printed = _seq_file_write_seq_plain(sf, seq);
      break;
    default:
      fprintf(stderr, "seq_file.c: called seq_file_write_seq() with invalid "
                      "file type (%s) [path: %s]\n", seq_file_get_type_str(sf),
                      sf->path);
      return 0;
  }

  sf->total_bases_passed += str_len;
  sf->write_state = WS_SEQ;

  return num_bytes_printed;
}

// Print quality
// Only FASTQ file types are allowed to call this function
size_t seq_file_write_qual(SeqFile *sf, const char *qual)
{
  if(sf->file_type != SEQ_FASTQ)
  {
    fprintf(stderr, "seq_file.c: called seq_file_write_qual() with invalid "
                    "file type (%s) [path: %s]\n", seq_file_get_type_str(sf),
                    sf->path);
    return 0;
  }

  if(sf->write_state == WS_BEGIN || sf->write_state == WS_NAME)
  {
    fprintf(stderr, "seq_file.c: writing in the wrong order (qual) [path: %s]\n",
            sf->path);
    exit(EXIT_FAILURE);
  }

  size_t str_len = strlen(qual);

  if(str_len == 0)
    return 0;

  size_t num_bytes_printed = 0;

  if(sf->write_state == WS_SEQ)
  {
    num_bytes_printed += seq_puts(sf, "\n+\n");
  }

  num_bytes_printed += _write(sf, qual, str_len);
  sf->write_state = WS_QUAL;

  return num_bytes_printed;
}
