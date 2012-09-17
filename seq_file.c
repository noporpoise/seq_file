/*
 seq_file.c
 project: seq_file
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
#include <ctype.h> // tolower
#include <zlib.h>

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
#define seq_puts(f,str) (size_t) \
((f)->plain_file != NULL ? (size_t)fputs((str), (f)->plain_file) \
                         : (size_t)gzputs((f)->gz_file, (str)))

// wrapper for fwrite/gzwrite
#define seq_write(f,str,len) (size_t) \
((f)->plain_file != NULL \
  ? (size_t)fwrite((str), sizeof(char), (size_t)(len), (f)->plain_file) \
  : (size_t)gzwrite((f)->gz_file, (str), (unsigned int)(len)))

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

  char fastq_ascii_offset; // defaults to 33

  enum SeqFileType file_type;

  // have we seen a '>' at the start of a line in a fasta file?
  // or a '@' in a fastq?
  // For 'plain' format files this is used to store the first char per entry
  char read_line_start;

  // name, index and bases-read/offset of current entry
  StrBuf *entry_name;
  unsigned long entry_index;
  
  unsigned long entry_offset, entry_offset_qual;

  // Whether an entry has been read in
  char entry_read, entry_read_qual;

  // Buffer for reading in bases in FASTQ files
  StrBuf *bases_buff;

  // Total bases read/written - initially 0
  unsigned long total_bases_passed;
  // Total bases skipped (not read through API) in file so far
  unsigned long total_bases_skipped;

  unsigned long line_number;

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

void _str_to_lower(char *str)
{
  for(; *str != '\0'; str++)
    *str = (char)tolower(*str);
}

char _contains_extension(const char *path, const char *ext, const size_t len)
{
  const char *tmp;

  for(tmp = path; (tmp = strstr(tmp, ext)) != NULL; tmp++)
  {
    // Check extension is followed by end-of-string or separator . or _
    if(*(tmp+len) == '\0' || *(tmp+len) == '.' || *(tmp+len) == '_')
    {
      return 1;
    }
  }

  return 0;
}

void seq_guess_filetype_from_path(const char *path, SeqFileType *file_type,
                                  char *zipped)
{
  #define ARRLEN 20

  const char* exts[ARRLEN] = {".fa",".fasta",
                              ".fq",".fastq",
                              ".faz",".fagz",".fa.gz",".fa.gzip",".fasta.gzip",
                              ".fqz",".fqgz",".fq.gz",".fq.gzip",".fastq.gzip",
                              ".txt",".txtgz",".txt.gz",".txt.gzip",
                              ".sam", ".bam"};

  const SeqFileType types[ARRLEN]
    = {SEQ_FASTA, SEQ_FASTA,
       SEQ_FASTQ, SEQ_FASTQ,
       SEQ_FASTA, SEQ_FASTA, SEQ_FASTA, SEQ_FASTA, SEQ_FASTA,
       SEQ_FASTQ, SEQ_FASTQ, SEQ_FASTQ, SEQ_FASTQ, SEQ_FASTQ,
       SEQ_PLAIN, SEQ_PLAIN, SEQ_PLAIN, SEQ_PLAIN,
       SEQ_SAM, SEQ_BAM};

  const char zips[ARRLEN] = {0, 0,
                             0, 0,
                             1, 1, 1, 1, 1,
                             1, 1, 1, 1, 1,
                             0, 1, 1, 1,
                             0, 0};

  int i;
  size_t strlens[ARRLEN];

  for(i = 0; i < ARRLEN; i++)
    strlens[i] = strlen(exts[i]);

  *file_type = SEQ_UNKNOWN;

  size_t path_len = strlen(path);

  char* strcpy = strdup(path);
  _str_to_lower(strcpy);

  // Check for an extension at the end
  for(i = 0; i < ARRLEN; i++)
  {
    if(path_len >= strlens[i] &&
       strcasecmp(path + path_len - strlens[i], exts[i]) == 0)
    {
      *file_type = types[i];
      *zipped = zips[i];
      free(strcpy);
      return;
    }
  }

  // Check for an extension anywhere in the filename
  for(i = 0; i < ARRLEN; i++)
  {
    if(_contains_extension(path, exts[i], strlens[i]))
    {
      *file_type = types[i];
      *zipped = zips[i];
      free(strcpy);
      return;
    }
  }
}

// Determines file type and opens necessary streams + mallocs memory
void _set_seq_filetype(SeqFile *sf)
{
  // Guess filetype from path
  SeqFileType file_type;
  char zipped;

  seq_guess_filetype_from_path(sf->path, &file_type, &zipped);

  if(file_type == SEQ_SAM)
  {
    // SAM
    sf->sam_file = samopen(sf->path, "r", 0);
    sf->file_type = SEQ_SAM;
    sf->bam = bam_init1();
    return;
  }
  else if(file_type == SEQ_BAM)
  {
    // BAM
    sf->sam_file = samopen(sf->path, "rb", 0);
    sf->file_type = SEQ_BAM;
    sf->bam = bam_init1();
    return;
  }

  // If not SAM or BAM, we can open it and determine its contents -
  // more reliable

  // Open file for the first time
  if(strcmp(sf->path, "-") == 0)
  {
    sf->gz_file = gzdopen(fileno(stdin), "r");
  }
  else
  {
    sf->gz_file = gzopen(sf->path, "r");
  }

  if(sf->gz_file == NULL)
  {
    return;
  }

  int first_char;

  // Move sf->line_number from 0 to 1 on first character
  // Then for each newline, line_number++

  do
  {
    first_char = gzgetc(sf->gz_file);
    sf->line_number++;
  } while (first_char != -1 && (first_char == '\n' || first_char == '\r'));

  if(first_char == -1)
  {
    fprintf(stderr, "seq_file.c warning: empty sequence file\n");
    return;
  }
  else if(first_char == '>')
  {
    // Reading FASTA
    sf->file_type = SEQ_FASTA;
    sf->read_line_start = 1;
  }
  else if(first_char == '@')
  {
    // Reading FASTQ
    sf->file_type = SEQ_FASTQ;
    sf->read_line_start = 1;
    sf->bases_buff = strbuf_new();
  }
  else if(is_base_char(first_char))
  {
    // Plain file
    sf->file_type = SEQ_PLAIN;
    sf->read_line_start = 0;

    if(gzungetc(first_char, sf->gz_file) == -1)
    {
      fprintf(stderr, "seq_file.c error: gzungetc failed\n");
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    fprintf(stderr, "seq_file.c warning: unknown filetype starting '%c'\n",
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

  sf->fastq_ascii_offset = 33;

  sf->file_type = SEQ_UNKNOWN;
  sf->read_line_start = 0;

  sf->entry_name = strbuf_new();
  sf->entry_index = 0;

  sf->entry_offset = 0;
  sf->entry_offset_qual = 0;

  sf->entry_read = 0;
  sf->entry_read_qual = 0;

  sf->bases_buff = NULL;

  sf->total_bases_passed = 0;
  sf->total_bases_skipped = 0;

  sf->line_number = 0;

  // For writing
  sf->plain_file = NULL;
  sf->line_wrap = 0;
  sf->curr_line_length = 0;
  sf->write_state = WS_READ_ONLY;

  return sf;
}

SeqFile* seq_file_open(const char* file_path)
{
  SeqFile* sf = _create_default_seq_file(file_path);
  _set_seq_filetype(sf);

  if(sf->file_type == SEQ_UNKNOWN)
  {
    free(sf);
    return NULL;
  }
  else
  {
    return sf;
  }
}

SeqFile* seq_file_open_filetype(const char* file_path,
                                SeqFileType file_type)
{
  SeqFile* sf = _create_default_seq_file(file_path);
  sf->file_type = file_type;

  if(file_type == SEQ_FASTQ)
  {
    sf->bases_buff = strbuf_new();
  }

  switch(file_type)
  {
    case SEQ_SAM:
      sf->sam_file = samopen(sf->path, "r", 0);
      sf->bam = bam_init1();
      break;
    case SEQ_BAM:
      sf->sam_file = samopen(sf->path, "rb", 0);
      sf->bam = bam_init1();
      break;
    case SEQ_FASTA:
    case SEQ_FASTQ:
    case SEQ_PLAIN:
      if(strcmp(sf->path, "-") == 0)
      {
        sf->gz_file = gzdopen(fileno(stdin), "r");
      }
      else
      {
        sf->gz_file = gzopen(sf->path, "r");
      }
      break;
    default:
      fprintf(stderr, "seq_file.c warning: invalid SeqFileType in "
                      "function seq_file_open_filetype()\n");
      free(sf);
      return NULL;
  }

  return sf;
}

// Open to write
// file_type must be FASTA, FASTQ
// set gzip to != 0 to turn on gzipping output
// if line_wrap is != 0, sequence lines are wrapped
SeqFile* seq_file_open_write(const char* file_path, SeqFileType file_type,
                             char gzip, unsigned long line_wrap)
{
  if(file_type == SEQ_SAM || file_type == SEQ_BAM)
  {
    fprintf(stderr, "seq_file.c error: cannot write to a SAM or BAM file\n");
    return NULL;
  }
  else if(file_type == SEQ_UNKNOWN)
  {
    fprintf(stderr, "seq_file.c error: cannot open file type SEQ_UNKNOWN\n");
    return NULL;
  }
  else if(file_type == SEQ_PLAIN && line_wrap != 0)
  {
    fprintf(stderr, "seq_file.c warning: cannot set line wrap with 'plain' "
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
    sf->line_number++;
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
    strbuf_free(sf->bases_buff);
  }

  if(sf->plain_file != NULL)
  {
    fclose(sf->plain_file);
  }

  strbuf_free(sf->entry_name);

  free(sf);

  return num_bytes_printed;
}

SeqFileType seq_file_get_type(const SeqFile* sf)
{
  return sf->file_type;
}

const char* seq_file_get_type_str(const SeqFile* sf)
{
  return seq_file_types[sf->file_type];
}

const char* seq_file_type_str(SeqFileType file_type,
                              char zipped)
{
  return zipped ? seq_file_types_zipped[file_type] : seq_file_types[file_type];
}

// Get a pointer to the file path
const char* seq_get_path(const SeqFile* sf)
{
  return sf->path;
}

// Set FASTQ ASCII offset (also applies to SAM/BAM)
void seq_set_fastq_ascii_offset(SeqFile *sf, char fastq_ascii_offset)
{
  sf->fastq_ascii_offset = fastq_ascii_offset;
}

// Get FASTQ ASCII offset (also applies to SAM/BAM)
char seq_get_fastq_ascii_offset(const SeqFile *sf)
{
  return sf->fastq_ascii_offset;
}

// Get the number of bases read/written so far
unsigned long seq_total_bases_passed(const SeqFile *sf)
{
  return sf->total_bases_passed;
}

// Get the total bases skipped (not read through API) in file so far
unsigned long seq_total_bases_skipped(const SeqFile *sf)
{
  return sf->total_bases_skipped;
}

unsigned long seq_curr_line_number(const SeqFile *sf)
{
  return sf->line_number;
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
    strbuf_gzreadline(sf->entry_name, sf->gz_file);
    strbuf_chomp(sf->entry_name);

    sf->line_number++;
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
      while((c = gzgetc(sf->gz_file)) != -1 && c != '\n' && c != '\r')
      {
        sf->total_bases_skipped++;
      }

      if(c == -1)
        return 0;
      else
        sf->line_number++; // Must have read a new line

      // Read through end of line chars
      while((c = gzgetc(sf->gz_file)) != -1 && (c == '\n' || c == '\r'))
      {
        sf->line_number++;
      }

      if(c == -1)
      {
        return 0;
      }
      else if(c != '>')
      {
        sf->total_bases_skipped++;
      }
    }
    while(c != '>');

    // Read name
    strbuf_gzreadline(sf->entry_name, sf->gz_file);
    strbuf_chomp(sf->entry_name);
    sf->line_number++;

    return 1;
  }
}

void _seq_read_fastq_sequence(SeqFile *sf)
{
  strbuf_reset(sf->bases_buff);

  int c;

  while((c = gzgetc(sf->gz_file)) != -1 && c != '+')
  {
    if(c != '\r' && c != '\n')
    {
      strbuf_append_char(sf->bases_buff, (char)c);
      strbuf_gzreadline(sf->bases_buff, sf->gz_file);
      strbuf_chomp(sf->bases_buff);
    }

    sf->line_number++;
  }

  if(c == -1)
  {
    fprintf(stderr, "seq_file.c: missing + in FASTQ [file: %s]\n", sf->path);
  }

  // Read to end of separator line
  if(c != '\r' && c != '\n')
  {
    strbuf_gzskip_line(sf->gz_file);
  }
}

char _seq_next_read_fastq(SeqFile *sf)
{
  if(sf->read_line_start)
  {
    // Read name
    strbuf_gzreadline(sf->entry_name, sf->gz_file);
    strbuf_chomp(sf->entry_name);
    sf->line_number++;

    // Read whole sequence
    _seq_read_fastq_sequence(sf);

    sf->read_line_start = 0;
    return 1;
  }
  else
  {
    int c;

    // Count bases not read in
    sf->total_bases_skipped += (strbuf_len(sf->bases_buff) - sf->entry_offset);

    // Skip over remaining quality values
    while(sf->entry_offset_qual < strbuf_len(sf->bases_buff))
    {
      if((c = gzgetc(sf->gz_file)) == -1)
        return 0;

      if(c != '\r' && c != '\n')
        sf->entry_offset_qual++;
      else
        sf->line_number++;
    }

    // Skip newlines
    while((c = gzgetc(sf->gz_file)) != -1 && (c == '\n' || c == '\r'))
      sf->line_number++;

    if(c == -1)
      return 0;

    if(c != '@')
    {
      fprintf(stderr, "seq_file.c: FASTQ header does not begin with '@' [%c]\n",
              c);
      return 0;
    }

    // Read name
    strbuf_gzreadline(sf->entry_name, sf->gz_file);
    strbuf_chomp(sf->entry_name);
    sf->line_number++;

    // Read whole sequence
    _seq_read_fastq_sequence(sf);

    return 1;
  }
}

char _seq_next_read_bam(SeqFile *sf)
{
  if(sf->entry_read)
  {
    // Count skipped bases
    sf->total_bases_skipped += (unsigned long)sf->bam->core.l_qseq -
                               sf->entry_offset;
  }

  if(samread(sf->sam_file, sf->bam) < 0)
    return 0;

  // Get name
  strbuf_append_str(sf->entry_name, bam1_qname(sf->bam));

  return 1;
}

char _seq_next_read_plain(SeqFile *sf)
{
  if(sf->read_line_start)
  {
    sf->read_line_start = 0;

    int c;

    while((c = gzgetc(sf->gz_file)) != -1 && c != '\r' && c != '\n')
    {
      sf->total_bases_skipped++;
    }

    if(c == -1)
      return 0;
    else
      sf->line_number++;

    return 1;
  }
  else
  {
    // Check if we can read a base
    int c = gzgetc(sf->gz_file);

    if(c == -1)
    {
      return 0;
    }
    else if(c == '\n' || c == '\r')
    {
      sf->line_number++;
    }
    else
    {
      sf->read_line_start = (char)c;
    }

    return 1;
  }
}

// Returns 1 on success 0 if no more to read
char seq_next_read(SeqFile *sf)
{
  strbuf_reset(sf->entry_name);

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

unsigned long seq_get_bases_read(SeqFile *sf)
{
  return sf->entry_offset;
}

unsigned long seq_get_quals_read(SeqFile *sf)
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
  
  while((next = gzgetc(sf->gz_file)) != -1 && (next == '\n' || next == '\r'))
    sf->line_number++;

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
    *c = (char)next;
    return 1;
  }
}

char _seq_read_base_fastq(SeqFile *sf, char *c)
{
  if(sf->entry_offset < strbuf_len(sf->bases_buff))
  {
    *c = strbuf_get_char(sf->bases_buff, sf->entry_offset);
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
  unsigned long query_len = (unsigned long)sf->bam->core.l_qseq;

  if(sf->entry_offset >= query_len)
    return 0;

  uint8_t *seq = bam1_seq(sf->bam);

  if(is_reversed)
  {
    unsigned long index = query_len - sf->entry_offset - 1;
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
  if(sf->read_line_start != 0)
  {
    *c = sf->read_line_start;
    sf->read_line_start = 0;
    return 1;
  }

  int next = gzgetc(sf->gz_file);

  if(next == -1)
  {
    return 0;
  }
  else if(next == '\r' || next == '\n')
  {
    sf->line_number++;
    return 0;
  }
  else
  {
    *c = (char)next;
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
  if(sf->entry_offset_qual >= strbuf_len(sf->bases_buff))
    return 0;

  int next;

  while((next = gzgetc(sf->gz_file)) != -1 && (next == '\n' || next == '\r'))
    sf->line_number++;

  if(next == -1)
  {
    fprintf(stderr, "seq_file.c: fastq file ended without finishing quality "
                    "scores [file: %s]\n", sf->path);
    return 0;
  }

  *c = (char)next;
  return 1;
}

char _seq_read_qual_bam(SeqFile *sf, char *c)
{
  char is_reversed = sf->bam->core.flag & 16;
  unsigned long query_len = (unsigned long)sf->bam->core.l_qseq;

  if(sf->entry_offset_qual >= query_len)
    return 0;

  uint8_t *seq = bam1_qual(sf->bam);
  unsigned long index;

  if(is_reversed)
    index = query_len - sf->entry_offset_qual - 1;
  else
    index = sf->entry_offset_qual;

  *c = sf->fastq_ascii_offset + seq[index];

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

// str must be at least k+1 bytes long.  Null-terminates str at position k+1.
// returns 1 on success, 0 otherwise
char seq_read_k_bases(SeqFile *sf, char* str, int k)
{
  // Check if we have read anything in
  if(!sf->entry_read)
    return 0;

  int i;

  for(i = 0; i < k; i++)
  {
    if(!seq_read_base(sf, str+i))
    {
      sf->entry_read = 0;
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

  int i;

  for(i = 0; i < k; i++)
  {
    if(!seq_read_qual(sf, str+i))
    {
      sf->entry_read_qual = 0;
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

char _seq_read_all_bases_bam(SeqFile *sf, StrBuf *sbuf)
{
  unsigned long qlen = (unsigned long)sf->bam->core.l_qseq;

  if(sf->entry_offset >= qlen)
    return 0;

  // Get reverse
  char is_reversed = sf->bam->core.flag & 16;

  strbuf_ensure_capacity(sbuf, qlen - sf->entry_offset);

  uint8_t *seq = bam1_seq(sf->bam);

  // read in and reverse complement (if needed)
  unsigned long i;
  for(i = sf->entry_offset; i < qlen; i++)
  {
    unsigned long index = (is_reversed ? i : qlen - i - 1);
    int8_t b = bam1_seqi(seq, index);
    char c = bam_nt16_rev_table[is_reversed ? seq_comp_table[b] : b];
    strbuf_append_char(sbuf, c);
  }

  return 1;
}

char _seq_read_all_bases_fasta(SeqFile *sf, StrBuf *sbuf)
{
  int c;

  while((c = gzgetc(sf->gz_file)) != -1 && c != '>')
  {
    if(c != '\r' && c != '\n')
    {
      strbuf_append_char(sbuf, (char)c);
      strbuf_gzreadline(sbuf, sf->gz_file);
      strbuf_chomp(sbuf);
    }

    sf->line_number++;
  }

  if(c == '>')
    sf->read_line_start = 1;

  return 1;
}

char _seq_read_all_bases_fastq(SeqFile *sf, StrBuf *sbuf)
{
  // Copy from buffer
  t_buf_pos len = sf->bases_buff->len - sf->entry_offset;
  strbuf_copy(sbuf, 0, sf->bases_buff, sf->entry_offset, len);

  return 1;
}

char _seq_read_all_bases_plain(SeqFile *sf, StrBuf *sbuf)
{
  t_buf_pos len = 0;

  if(sf->read_line_start)
  {
    strbuf_append_char(sbuf, sf->read_line_start);
    sf->read_line_start = 0;
    len++;
  }

  len += strbuf_gzreadline(sbuf, sf->gz_file);
  strbuf_chomp(sbuf);
  sf->line_number++;

  return (len > 0);
}

// returns 1 on success, 0 otherwise
char seq_read_all_bases(SeqFile *sf, StrBuf *sbuf)
{
  // Check if we have read anything in
  if(!sf->entry_read)
    return 0;

  strbuf_reset(sbuf);

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

  sf->entry_offset += strbuf_len(sbuf);
  sf->total_bases_passed += strbuf_len(sbuf);

  // read has been exhausted
  sf->entry_read = 0;

  return success;
}


/*
 Read the rest of a read's quality scores
*/

char _seq_read_all_quals_bam(SeqFile *sf, StrBuf *sbuf)
{
  unsigned long qlen = (unsigned long)sf->bam->core.l_qseq;

  if(sf->entry_offset_qual >= qlen)
    return 0;

  // Get reverse
  char is_reversed = sf->bam->core.flag & 16;

  strbuf_ensure_capacity(sbuf, qlen - sf->entry_offset);

  uint8_t *seq = bam1_qual(sf->bam);

  // read in and reverse complement (if needed)
  unsigned long i;
  for(i = sf->entry_offset; i < qlen; i++)
  {
    char c = sf->fastq_ascii_offset + seq[is_reversed ? i : qlen - i - 1];
    strbuf_append_char(sbuf, c);
  }

  return 1;
}

char _seq_read_all_quals_fastq(SeqFile *sf, StrBuf *sbuf)
{
  if(sf->entry_offset_qual >= strbuf_len(sf->bases_buff))
    return 0;

  // Expect the same number of quality scores as bases
  t_buf_pos expected_len = strbuf_len(sf->bases_buff) -
                           sf->entry_offset_qual;

  int next = -1;
  t_buf_pos i;

  for(i = 0; i < expected_len && (next = gzgetc(sf->gz_file)) != -1; i++)
  {
    if(next != '\r' && next != '\n')
    {
      strbuf_append_char(sbuf, (char)next);
    }
    else
      sf->line_number++;
  }

  if(next == -1)
  {
    fprintf(stderr, "seq_file.c: fastq file ended without finishing quality "
                    "scores (FASTQ) [file: %s]\n", sf->path);
  }

  return 1;
}

// returns 1 on success, 0 otherwise
char seq_read_all_quals(SeqFile *sf, StrBuf *sbuf)
{
  // Check if we have read anything in
  if(!sf->entry_read_qual)
    return 0;

  strbuf_reset(sbuf);

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
  sf->entry_offset_qual += strbuf_len(sbuf);
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
  else
  {
    num_bytes_printed += seq_puts(sf, "\n>");
    sf->line_number++;
  }

  num_bytes_printed += seq_puts(sf, name);

  return num_bytes_printed;
}

unsigned long _seq_file_write_name_fastq(SeqFile *sf, const char *name)
{
  unsigned long num_bytes_printed = 0;

  if(sf->write_state == WS_BEGIN)
  {
    num_bytes_printed += seq_puts(sf, "@");
  }
  else if(sf->write_state == WS_NAME)
  {
    num_bytes_printed += seq_puts(sf, "\n\n+\n\n@");
  }
  else if(sf->write_state == WS_QUAL)
  {
    num_bytes_printed += seq_puts(sf, "\n@");
    sf->line_number++;
  }
  else if(sf->write_state == WS_SEQ)
  {
    fprintf(stderr, "seq_file.c: writing in the wrong order (name) [path: %s]\n",
            sf->path);
    exit(EXIT_FAILURE);
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

size_t _write_wrapped(SeqFile *sf, const char *str, size_t str_len)
{
  size_t num_bytes_printed = 0;

  if(sf->curr_line_length == sf->line_wrap)
  {
    sf->curr_line_length = 0;
    num_bytes_printed += seq_puts(sf, "\n");
    sf->line_number++;
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
  sf->line_number++;

  size_t offset;

  for(offset = bytes_to_print; offset < str_len; offset += sf->line_wrap)
  {
    bytes_to_print = MIN(str_len - offset, sf->line_wrap);
    num_bytes_printed += (size_t)seq_write(sf, str + offset, bytes_to_print);

    if(bytes_to_print < sf->line_wrap)
    {
      num_bytes_printed += seq_puts(sf, "\n");
      sf->line_number++;
    }
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
    sf->line_number++;
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
    sf->line_number++;
  }

  num_bytes_printed += _write(sf, seq, str_len);

  return num_bytes_printed;
}

size_t _seq_file_write_seq_plain(SeqFile *sf, const char *seq)
{
  size_t num_bytes_printed = 0;

  if(sf->write_state != WS_BEGIN)
  {
    num_bytes_printed += seq_puts(sf, "\n");
    sf->line_number++;
  }

  num_bytes_printed += seq_puts(sf, seq);

  return num_bytes_printed;
}

size_t seq_file_write_seq(SeqFile *sf, const char *seq)
{
  size_t str_len = strlen(seq);

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
    sf->line_number += 2;
  }

  num_bytes_printed += _write(sf, qual, str_len);
  sf->write_state = WS_QUAL;

  return num_bytes_printed;
}
