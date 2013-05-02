
#ifndef _SEQ_FILE_HEADER
#define _SEQ_FILE_HEADER

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h> // strcasecmp
#include <ctype.h>
#include <limits.h>
#include <zlib.h>

#include "hts.h"
#include "sam.h"
#include "stream_buffer.h"

typedef enum
{
  IS_ERROR, IS_UNKNOWN,
  IS_SEQ_UNKNOWN, // could be FASTQ, FASTQ or plain
  IS_SAM, IS_BAM,
  IS_FASTQ, IS_FASTA, IS_PLAIN
} seqtype_t;

typedef struct seq_file_t seq_file_t;
typedef struct read_t read_t;

struct seq_file_t
{
  char *path;
  FILE *f_file;
  gzFile gz_file;
  samFile *s_file;
  bam_hdr_t *bam_header;
  int (*read)(seq_file_t *sf, read_t *r);
  buffer_t in;
  int nextc;
  seqtype_t format;
};

struct read_t
{
  buffer_t name, seq, qual;
  bam1_t *bam;
  char from_sam; // from sam or bam
};

#define seq_is_bam(sf) ((sf)->format == IS_BAM)
#define seq_is_sam(sf) ((sf)->format == IS_SAM)
#define seq_use_gzip(sf) ((sf)->gz_file != NULL)

// The following require a read to have been read successfully first
// using seq_read
#define seq_is_fastq(sf) ((sf)->format == IS_FASTQ)
#define seq_is_fasta(sf) ((sf)->format == IS_FASTA)
#define seq_is_plain(sf) ((sf)->format == IS_PLAIN)

#define seq_get_path(sf) ((sf)->path)

// return 1 on success, 0 on eof, -1 if partially read / syntax error
#define seq_read(sf,r) (sf)->read(sf,r)

// File format information (http://en.wikipedia.org/wiki/FASTQ_format)
static const char * const FASTQ_FORMATS[]
  = {"Sanger / Illumina 1.9+ (Phred+33)", // range: [0,71] "catch all / unknown"
     "Sanger (Phred+33)", // range: [0,40]
     "Solexa (Solexa+64)", // range: [-5,40]
     "Illumina 1.3+ (Phred+64)", // range: [0,40]
     "Illumina 1.5+ (Phred+64)", // range: [3,40]
     "Illumina 1.8+ (Phred+33)"}; // range: [0,41]

static const int FASTQ_MIN[6]    = { 33, 33, 59, 64, 67, 33};
static const int FASTQ_MAX[6]    = {126, 73,104,104,104, 74};
static const int FASTQ_OFFSET[6] = { 33, 33, 64, 64, 64, 33};

//

// file could be sam,bam,FASTA,FASTQ,txt (+gzip)

#define DEFAULT_BUFSIZE (1<<20)

static const int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 9, 14,
                                           1, 6, 5, 13, 3, 11, 7, 15 };


// Read a sam/bam file
static inline int sread_s(seq_file_t *sf, read_t *read)
{
  read->name.end = read->seq.end = read->qual.end = 0;
  read->name.b[0] = read->seq.b[0] = read->qual.b[0] = '\0';

  if(sam_read1(sf->s_file, sf->bam_header, read->bam) < 0) return 0;

  char *str = bam_get_qname(read->bam);
  buffer_append_str(&read->name, str);

  size_t qlen = read->bam->core.l_qseq;
  buffer_ensure_capacity(&read->seq, qlen);
  buffer_ensure_capacity(&read->qual, qlen);
  uint8_t *bamseq = bam_get_seq(read->bam);
  uint8_t *bamqual = bam_get_qual(read->bam);

  size_t i, j;
  if(bam_is_rev(read->bam))
  {
    for(i = 0, j = qlen - 1; i < qlen; i++, j--)
    {
      int8_t b = bam_seqi(bamseq, j);
      read->seq.b[i] = seq_nt16_str[seq_comp_table[b]];
      read->qual.b[i] = 33 + bamqual[j];
    }
  }
  else
  {
    for(i = 0; i < qlen; i++)
    {
      int8_t b = bam_seqi(bamseq, i);
      read->seq.b[i] = seq_nt16_str[b];
      read->qual.b[i] = 33 + bamqual[i];
    }
  }

  read->seq.end = read->qual.end = qlen;
  read->seq.b[qlen] = read->qual.b[qlen] = 0;
  read->from_sam = 1;

  return 1;
}


#define RESET_READING(sf,read) do { \
  read->name.end = read->seq.end = read->qual.end = 0;        \
  read->name.b[0] = read->seq.b[0] = read->qual.b[0] = '\0';  \
  read->from_sam = 0;                                         \
  sf->nextc = -1;                                             \
} while(0)

#define _func_read_fastq(_read_fastq,__getc,__readline) \
  static inline int _read_fastq(seq_file_t *sf, read_t *read)                  \
  {                                                                            \
    int c = (sf->nextc != -1 ? sf->nextc : __getc(sf));                        \
    RESET_READING(sf,read);                                                    \
                                                                               \
    if(c == -1) return 0;                                                      \
    if(c != '@' || __readline(sf, read->name) == 0) return -1;                 \
    buffer_chomp(&(read->name));                                               \
                                                                               \
    while((c = __getc(sf)) != '+') {                                           \
      if(c == -1) return -1;                                                   \
      if(c != '\r' && c != '\n') {                                             \
        buffer_append_char(&read->seq,c);                                      \
        if(__readline(sf, read->seq) == 0) return -1;                          \
        buffer_chomp(&(read->seq));                                            \
      }                                                                        \
    }                                                                          \
    while((c = __getc(sf)) != -1 && c != '\n');                                \
    if(c == -1) return -1;                                                     \
    do {                                                                       \
      if(__readline(sf,read->qual) > 0) buffer_chomp(&(read->qual));           \
      else return 1;                                                           \
    } while(read->qual.end < read->seq.end);                                   \
    while((c = __getc(sf)) != -1 && c != '@');                                 \
    sf->nextc = c;                                                             \
    return 1;                                                                  \
  }

#define _func_read_fasta(_read_fasta,__getc,__readline) \
  static inline int _read_fasta(seq_file_t *sf, read_t *read)                  \
  {                                                                            \
    int c = (sf->nextc != -1 ? sf->nextc : __getc(sf));                        \
    RESET_READING(sf,read);                                                    \
                                                                               \
    if(c == -1) return 0;                                                      \
    if(c != '>' || __readline(sf, read->name) == 0) return -1;                 \
    buffer_chomp(&(read->name));                                               \
                                                                               \
    while((c = __getc(sf)) != '>') {                                           \
      if(c == -1) return 1;                                                    \
      if(c != '\r' && c != '\n') {                                             \
        buffer_append_char(&read->seq,c);                                      \
        int r = __readline(sf, read->seq);                                     \
        buffer_chomp(&(read->seq));                                            \
        if(r <= 0) return 1;                                                   \
      }                                                                        \
    }                                                                          \
    sf->nextc = c;                                                             \
    return 1;                                                                  \
  }

#define _func_read_plain(_read_plain,__getc,__readline,__skipline)             \
  static inline int _read_plain(seq_file_t *sf, read_t *read)                  \
  {                                                                            \
    int c = (sf->nextc != -1 ? sf->nextc : __getc(sf));                        \
    RESET_READING(sf,read);                                                    \
    while(c != -1 && isspace(c) && c != '\n') { __skipline(sf); }              \
    if(c == -1) return 0;                                                      \
    buffer_append_char(&read->seq, c);                                         \
    __readline(sf, read->seq);                                                 \
    buffer_chomp(&(read->seq));                                                \
    return 1;                                                                  \
  }

#define _func_read_unknown(_read_unknown,__getc,__skipline,_fastq,_fasta,_plain)\
  static inline int _read_unknown(seq_file_t *sf, read_t *read)                \
  {                                                                            \
    RESET_READING(sf, read);                                                   \
    int c;                                                                     \
    while((c = __getc(sf)) != -1 && isspace(c)) {                              \
      if(c != '\n') { __skipline(sf); }                                        \
    }                                                                          \
    if(c == -1) return 0;                                                      \
    if(c == '@') { sf->format = IS_FASTQ; sf->read = _fastq; }                 \
    else if(c == '>') { sf->format = IS_FASTA; sf->read = _fasta; }            \
    else { sf->format = IS_PLAIN; sf->read = _plain; }                         \
    sf->nextc = c;                                                             \
    return sf->read(sf,read);                                                  \
  }

// perform reading on seq_file_t

// getc on seq_file_t
#define _sf_gzgetc(sf)              gzgetc(sf->gz_file)
#define _sf_gzgetc_buf(sf)          gzgetc_buf(sf->gz_file,&sf->in)
#define _sf_fgetc(sf)               fgetc(sf->f_file)
#define _sf_fgetc_buf(sf)           fgetc_buf(sf->f_file,&sf->in)

// readline on seq_file_t
#define _sf_gzreadline(sf,buf)      gzreadline(sf->gz_file,&buf.b,&buf.end,&buf.size)
#define _sf_gzreadline_buf(sf,buf)  gzreadline_buf(sf->gz_file,&sf->in,&buf.b,&buf.end,&buf.size)
#define _sf_freadline(sf,buf)       freadline(sf->f_file,&buf.b,&buf.end,&buf.size)
#define _sf_freadline_buf(sf,buf)   freadline_buf(sf->f_file,&sf->in,&buf.b,&buf.end,&buf.size)

// skipline on seq_file_t
#define _sf_gzskipline(sf)          gzskipline(sf->gz_file)
#define _sf_gzskipline_buf(sf)      gzskipline_buf(sf->gz_file,&sf->in)
#define _sf_fskipline(sf)           fskipline(sf->f_file)
#define _sf_fskipline_buf(sf)       fskipline_buf(sf->f_file,&sf->in)

// Read FASTQ
_func_read_fastq(_seq_read_fastq_f,      _sf_fgetc,      _sf_freadline)
_func_read_fastq(_seq_read_fastq_gz,     _sf_gzgetc,     _sf_gzreadline)
_func_read_fastq(_seq_read_fastq_f_buf,  _sf_fgetc_buf,  _sf_freadline_buf)
_func_read_fastq(_seq_read_fastq_gz_buf, _sf_gzgetc_buf, _sf_gzreadline_buf)

// Read FASTA
_func_read_fasta(_seq_read_fasta_f,      _sf_fgetc,      _sf_freadline)
_func_read_fasta(_seq_read_fasta_gz,     _sf_gzgetc,     _sf_gzreadline)
_func_read_fasta(_seq_read_fasta_f_buf,  _sf_fgetc_buf,  _sf_freadline_buf)
_func_read_fasta(_seq_read_fasta_gz_buf, _sf_gzgetc_buf, _sf_gzreadline_buf)

// Read plain
_func_read_plain(_seq_read_plain_f,      _sf_fgetc,      _sf_freadline,      _sf_fskipline)
_func_read_plain(_seq_read_plain_gz,     _sf_gzgetc,     _sf_gzreadline,     _sf_gzskipline)
_func_read_plain(_seq_read_plain_f_buf,  _sf_fgetc_buf,  _sf_freadline_buf,  _sf_fskipline_buf)
_func_read_plain(_seq_read_plain_gz_buf, _sf_gzgetc_buf, _sf_gzreadline_buf, _sf_gzskipline_buf)

// Read first entry
_func_read_unknown(_seq_read_unknown_f,      _sf_fgetc,  _sf_fskipline,  _seq_read_fastq_f,      _seq_read_fasta_f,      _seq_read_plain_f)
_func_read_unknown(_seq_read_unknown_gz,     _sf_gzgetc, _sf_gzskipline, _seq_read_fastq_gz,     _seq_read_fasta_gz,     _seq_read_plain_gz)
_func_read_unknown(_seq_read_unknown_f_buf,  _sf_fgetc,  _sf_fskipline,  _seq_read_fastq_f_buf,  _seq_read_fasta_f_buf,  _seq_read_plain_f_buf)
_func_read_unknown(_seq_read_unknown_gz_buf, _sf_gzgetc, _sf_gzskipline, _seq_read_fastq_gz_buf, _seq_read_fasta_gz_buf, _seq_read_plain_gz_buf)

// Create and destroy read structs
static inline read_t* seq_read_alloc()
{
  read_t *r = calloc(1, sizeof(read_t));
  buffer_init(&r->name, 512);
  buffer_init(&r->seq, DEFAULT_BUFSIZE);
  buffer_init(&r->qual, DEFAULT_BUFSIZE);
  r->bam = bam_init1();
  r->from_sam = 0;
  return r;
}

static inline void seq_read_destroy(read_t *r)
{
  free(r->name.b);
  free(r->seq.b);
  free(r->qual.b);
  free(r->bam);
  free(r);
}

#define seq_file_init(sf) do { \
    sf->gz_file = NULL;                                                        \
    sf->f_file = NULL;                                                         \
    sf->s_file = NULL;                                                         \
    sf->nextc = -1;                                                            \
    sf->in.size = sf->in.begin = sf->in.end = 0;                               \
    sf->path = sf->in.b = NULL;                                                \
    sf->format = IS_UNKNOWN;                                                   \
  } while(0)

// I have removed gzbuffer for now since it causes linking errors on systems
// that are not set up properly.  Feel free to uncomment and remove empty
// definition
#if defined(ZLIB_VERNUM) && ZLIB_VERNUM >= 0x1240
//#define SET_ZLIB_BUFFER(sf,s) gzbuffer((sf)->gz_file, DEFAULT_BUFSIZE)
#define SET_ZLIB_BUFFER(sf,s)
#else
#define SET_ZLIB_BUFFER(sf,s)
#endif

#define seq_setup(sf,use_zlib,buf_size) do {\
  sf->in.size = buf_size;                                                      \
  if(sf->in.size > 0 && (sf->in.b = malloc(sf->in.size)) == NULL) {            \
    free(sf); return NULL;                                                     \
  }                                                                            \
  if((sf)->in.size == 0 && use_zlib) { SET_ZLIB_BUFFER(sf,DEFAULT_BUFSIZE); }  \
  if(use_zlib) sf->read = sf->in.size > 0 ? _seq_read_unknown_gz_buf : _seq_read_unknown_gz; \
  else         sf->read = sf->in.size > 0 ? _seq_read_unknown_f_buf  : _seq_read_unknown_f;  \
} while(0)

// Guess file type from file path or contents
#define EXT_ARRLEN 28

// Reports 
static inline seqtype_t _guess_filetype_from_filename(const char *path)
{
  size_t plen = strlen(path);
  const char* exts[EXT_ARRLEN]
    = {".fa", ".fasta", ".fsa", ".fsa.gz", "fsa.gzip", // FASTA
       ".faz", ".fagz", ".fa.gz", ".fa.gzip", ".fastaz", ".fasta.gzip",
       ".fq", ".fastq", ".fsq", ".fsq.gz", "fsq.gzip", // FASTQ
       ".fqz", ".fqgz", ".fq.gz", ".fq.gzip", ".fastqz", ".fastq.gzip",
       ".txt", ".txtgz", ".txt.gz", ".txt.gzip", // Plain
       ".sam", ".bam"}; // SAM / BAM

  const seqtype_t types[EXT_ARRLEN]
    = {IS_FASTA, IS_FASTA, IS_FASTA, IS_FASTA, IS_FASTA,
       IS_FASTA, IS_FASTA, IS_FASTA, IS_FASTA,
       IS_FASTQ, IS_FASTQ, IS_FASTQ, IS_FASTQ, IS_FASTQ,
       IS_FASTQ, IS_FASTQ, IS_FASTQ, IS_FASTQ,
       IS_PLAIN, IS_PLAIN, IS_PLAIN, IS_PLAIN,
       IS_SAM, IS_BAM};

  size_t extlens[EXT_ARRLEN];
  size_t i;
  for(i = 0; i < EXT_ARRLEN; i++)
    extlens[i] = strlen(exts[i]);

  for(i = 0; i < EXT_ARRLEN; i++)
    if(extlens[i] <= plen && strcasecmp(path+plen-extlens[i], exts[i]) == 0)
      return types[i];

  return IS_UNKNOWN;
}

// str should point to a tab character
static inline char _is_num_column(char *str, char ** ptr)
{
  int i;
  for(i = 1; str[i] >= '0' && str[i] <= '9'; i++);
  *ptr = str+i;
  return (i > 1 && str[i] == '\t');
}

// str should point to a tab character
static inline char _is_cigar_column(char *str, char ** ptr)
{
  int i = 1;
  if(str[i] == '*' && str[i+1] == '\t') {
    *ptr = str + 2;
    return 1;
  }
  while(str[i] != '\0') {
    if(str[i] >= '0' && str[i] <= '9') i++;
    else return 0;
    while(str[i] >= '0' && str[i] <= '9') i++;
    if(strchr("MIDNSHP=X", str[i])) i++;
    else break;
  }
  *ptr = str+i;
  return (i > 1 && str[i] == '\t');
}

#define ISUPPER(x) ((x) >= 'A' && (x) <= 'Z')

static inline seqtype_t _guess_filetype_from_content(const char *buf, int len)
{
  // printf("buf: '%s'\n", buf);
  if(len >= 4 && strncmp(buf, "@HD\t", 4) == 0) return IS_SAM;
  if(len >= 3 && strncmp(buf, "BAM", 3) == 0) return IS_BAM;

  // Test for @UU\tUU: where U is any uppercase char
  if(buf[0] == '@' && ISUPPER(buf[1]) && ISUPPER(buf[2]) &&
     buf[3] == '\t' && ISUPPER(buf[4]) && ISUPPER(buf[5])) return IS_SAM;

  // Test for a sam entry with no headers
  char *tmp = strchr(buf,'\t');
  if(tmp != NULL) {
    if(!_is_num_column(tmp, &tmp)) return IS_SEQ_UNKNOWN;
    tmp = strchr(tmp+1, '\t');
    if(!_is_num_column(tmp, &tmp)) return IS_SEQ_UNKNOWN;
    if(!_is_num_column(tmp, &tmp)) return IS_SEQ_UNKNOWN;
    if(!_is_cigar_column(tmp, &tmp)) return IS_SEQ_UNKNOWN;
    return IS_SAM;
  }
  return IS_SEQ_UNKNOWN;
}

#define GUESS_BUFLEN 1000

static inline seqtype_t _guess_filetype_from_content_using_path(const char *path)
{
  char buf[GUESS_BUFLEN];
  gzFile gz = gzopen(path, "r");
  if(gz == NULL) return IS_ERROR;
  int len = gzread(gz, buf, GUESS_BUFLEN);
  gzclose(gz);

  return _guess_filetype_from_content(buf, len);
}

static inline seqtype_t _guess_filetype_from_content_using_fh(FILE *fh)
{
  char buf[GUESS_BUFLEN];
  int c, i, len = 0;

  // Read first line
  while(len < GUESS_BUFLEN-1 && (c = fgetc(fh)) != -1)
  {
    buf[len++] = c;
    if(c == '\n' || c == '\r') break;
  }
  buf[len] = '\0';

  // Unget sequence
  for(i = len-1; i >= 0 && ungetc(buf[i], fh) != EOF; i--);

  // Check if we were able to restore the buffer to the FILE stream
  if(i >= 0)
  {
    fprintf(stderr, "%s:%c:Error: Cannot push back onto FILE stream\n",
            __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  return _guess_filetype_from_content(buf, len);
}

// format can be: IS_SAM, IS_BAM, IS_SEQ_UNKNOWN, IS_FASTA, IS_FASTQ, IS_PLAIN
// format cannot be IS_ERROR, IS_UNKNOWN
static inline seq_file_t* seq_open2(const char *p, seqtype_t format,
                                    char use_zlib, size_t buf_size)
{
  seq_file_t *sf = calloc(1, sizeof(seq_file_t));
  seq_file_init(sf);
  sf->path = strdup(p);
  sf->format = format;

  if(format == IS_ERROR || format == IS_UNKNOWN) {
    fprintf(stderr, "[%s:%i] Error: format invalid\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  if(format == IS_SAM || format == IS_BAM)
  {
    if((sf->s_file = sam_open(p, format == IS_SAM ? "rs" : "rb", 0)) == NULL) {
      free(sf->path);
      free(sf);
      return NULL;
    }
    sf->bam_header = sam_hdr_read(sf->s_file);
    sf->read = sread_s;
    return sf;
  }

  if(( use_zlib && ((sf->gz_file = gzopen(p, "r")) == NULL)) ||
     (!use_zlib && ((sf->f_file  =  fopen(p, "r")) == NULL)))
  {
    free(sf->path);
    free(sf);
    return NULL;
  }

  seq_setup(sf, use_zlib, buf_size);

  return sf;
}

// Returns pointer to new seq_file_t on success, seq_close will close the fh,
// so you shouldn't call fclose(fh)
// Returns NULL on error, in which case FILE will not have been closed (caller
// should then call fclose(fh))
static inline seq_file_t* seq_open_fh2(FILE *fh, char sam_bam,
                                       char use_zlib, size_t buf_size)
{
  seq_file_t *sf = calloc(1, sizeof(seq_file_t));
  seq_file_init(sf);

  if(sam_bam == 1 || sam_bam == 2)
  {
    fprintf(stderr, "%s:%i:Error: opening SAM/BAM from a FILE stream not "
                    "implemented yet, sorry\n", __FILE__, __LINE__);
    free(sf);
    return NULL;

    /*
    if((sf->s_file = sam_open("-", sam_bam == 1 ? "rs" : "rb", 0)) == NULL) {
      free(sf);
      return NULL;
    }
    sf->bam_header = sam_hdr_read(sf->s_file);
    sf->read = sread_s;
    */
  }
  else
  {
    if(!use_zlib) sf->f_file = fh;
    else if((sf->gz_file = gzdopen(fileno(fh), "r")) == NULL) {
      free(sf);
      return NULL;
    }
    seq_setup(sf, use_zlib, buf_size);
  }

  sf->path = strdup("-");
  return sf;
}

static inline seq_file_t* seq_open_fh(FILE *fh, char buffered)
{
  seqtype_t type = _guess_filetype_from_content_using_fh(fh);

  // sam -> 1, bam -> 2, other -> 0
  char sam_bam = (type == IS_SAM ? 1 : (type == IS_BAM ? 2 : 0));

  char gzip = buffered ? 1 : 0;
  size_t bufsize = buffered ? DEFAULT_BUFSIZE : 0;
  return seq_open_fh2(fh, sam_bam, gzip, bufsize);
}

static inline seq_file_t* seq_open(const char *p)
{
  if(strcmp(p,"-") == 0) return seq_open_fh(stdin, 0);
  seqtype_t format = _guess_filetype_from_filename(p);

  if(format == IS_UNKNOWN) {
    format = _guess_filetype_from_content_using_path(p);
    if(format == IS_ERROR) return NULL;
  }

  // If not sam or bam, don't assume FASTA/FASTQ until we've seen the contents
  if(format != IS_SAM && format != IS_BAM) format = IS_SEQ_UNKNOWN;

  // Read as gzip file with default buffer size
  return seq_open2(p, format, 1, DEFAULT_BUFSIZE);
}

static inline void seq_close(seq_file_t *sf)
{
  if(sf->f_file != NULL) { fclose(sf->f_file); }
  else if(sf->gz_file != NULL) { gzclose(sf->gz_file); }
  else if(sf->s_file != NULL)
  {
    sam_close(sf->s_file);
    free(sf->bam_header);
  }
  if(sf->in.size != 0) { free(sf->in.b); }
  if(sf->path != NULL) { free(sf->path); }
}

// Get min and max quality values in the first `num` bases of a file.
// Returns 0 if no qual scores, 1 on success, -1 if read error
static inline int seq_get_qual_limits(const char *path, size_t num,
                                      int *minptr, int *maxptr)
{
  seq_file_t *sf = seq_open(path);
  if(sf == NULL) return -1;
  read_t *r = seq_read_alloc();
  int i, limit, min = INT_MAX, max = 0;
  size_t count = 0;

  while(count < num && seq_read(sf,r))
  {
    limit = (r->qual.end < num-count ? r->qual.end : num-count);
    for(i = 0; i < limit; i++)
    {
      char q = r->qual.b[i];
      if(q > max) max = q;
      if(q < min) min = q;
    }
    count += limit;
  }

  seq_close(sf);

  if(count > 0) {
    *minptr = min;
    *maxptr = max;
  }
  return (count > 0);
}

// Returns -1 on error
// Returns 0 if not in FASTQ format/ not recognisable (offset:33, min:33, max:104)
static inline int seq_guess_fastq_format(const char *path)
{
  // Detect fastq offset
  int min_qual = INT_MAX, max_qual = 0;

  // 1000 is the number of quality scores to read
  if(seq_get_qual_limits(path, 500, &min_qual, &max_qual) <= 0) {
    return -1;
  }

  // See: http://en.wikipedia.org/wiki/FASTQ_format
  // Usually expect 0,40, but new software can report 41, so using <= MAX+1
  if(min_qual >= 33 && max_qual <= 74) return 1; // sanger
  else if(min_qual >= 33 && max_qual <= 75) return 5; // Illumina 1.8+
  else if(min_qual >= 67 && max_qual <= 105) return 4; // Illumina 1.5+
  else if(min_qual >= 64 && max_qual <= 105) return 3; // Illumina 1.3+
  else if(min_qual >= 59 && max_qual <= 105) return 2; // Solexa
  else return 0; // Unknown, assume 33 offset max value 104
}

static inline char _seq_read_looks_valid(read_t *r, const char *alphabet)
{
  size_t i;
  if(r->qual.end != 0) {
    if(r->qual.end != r->seq.end) return 0;
    for(i = 0; i < r->seq.end; i++) {
      char b = tolower(r->seq.b[i]);
      char q = r->qual.b[i];
      if(strchr(alphabet, b) == NULL) return 0;
      if(q < 33 || q > 105) return 0;
    }
  }
  else {
    for(i = 0; i < r->seq.end; i++) {
      char b = tolower(r->seq.b[i]);
      if(strchr(alphabet, b) == NULL) return 0;
    }
  }
  return 1;
}

#define seq_read_looks_valid_dna(r) _seq_read_looks_valid(r,"acgtn")
#define seq_read_looks_valid_rna(r) _seq_read_looks_valid(r,"acgun")
#define seq_read_looks_valid_protein(r) _seq_read_looks_valid(r,"acdefghiklmnopqrstuvwy")

static inline char _seq_char_complement(char c) {
  switch(c) {
    case 'a': return 't'; case 'A': return 'T';
    case 'c': return 'g'; case 'C': return 'G';
    case 'g': return 'c'; case 'G': return 'C';
    case 't': return 'a'; case 'T': return 'A';
    case 'n': return 'n'; case 'N': return 'N';
    default: return c;
  }
}

static inline void seq_read_reverse_complement(read_t *r)
{
  size_t i, j;
  char swap;

  if(r->qual.end > 0)
  {
    // Force quality score length to match seq length
    if(r->qual.end < r->seq.end) {
      buffer_ensure_capacity(&(r->qual), r->seq.end);
      for(i = r->qual.end; i < r->seq.end; i++) r->qual.b[i] = '.';
    }
    r->qual.b[r->qual.end = r->seq.end] = '\0';
  }

  if(r->seq.end == 0) return;
  if(r->seq.end == 1){ r->seq.b[0] = _seq_char_complement(r->seq.b[0]); return; }

  for(i=0, j=r->seq.end-1; i <= j; i++, j--) {
    swap = r->seq.b[i];
    r->seq.b[i] = _seq_char_complement(r->seq.b[j]);
    r->seq.b[j] = _seq_char_complement(swap);
  }

  if(r->qual.end > 0)
  {
    for(i=0, j=r->qual.end-1; i <= j; i++, j--) {
      swap = r->qual.b[i];
      r->qual.b[i] = r->qual.b[j];
      r->qual.b[j] = swap;
    }
  }
}

#define _seq_print_wrap(fh,str,len,wrap,i,j,_putc) do { \
    for(i=0,j=0;i<len;i++,j++) { \
      if(j==wrap) { _putc((fh),'\n'); j = 0; } \
      _putc((fh),str[i]); \
    } \
  } while(0)

#define _seq_print_fasta(fname,ftype,_printf,_putc)                            \
  static inline void fname(const read_t *r, ftype fh, int linewrap) {          \
    if(linewrap == 0)  _printf(fh, ">%s\n%s\n", r->name.b, r->seq.b);          \
    else {                                                                     \
      size_t i; int j;                                                         \
      _printf(fh, ">%s\n", r->name.b);                                         \
      _seq_print_wrap(fh, r->seq.b, r->seq.end, linewrap, i, j, _putc);        \
      _putc(fh, '\n');                                                         \
    }                                                                          \
  }                                                                            \

_seq_print_fasta(seq_print_fasta,FILE*,fprintf,fputc2)
_seq_print_fasta(seq_gzprint_fasta,gzFile,gzprintf,gzputc2)

#define _seq_print_fastq(fname,ftype,_printf,_putc)                            \
    static inline void fname(const read_t *r, ftype fh, int linewrap) {        \
    size_t i, qlimit = (r->qual.end < r->seq.end ? r->qual.end : r->seq.end);  \
    int j;                                                                     \
    if(linewrap <= 0) {                                                        \
      _printf(fh, "@%s\n%s\n+\n%.*s",r->name.b,r->seq.b,(int)qlimit,r->qual.b);\
      for(i = r->qual.end; i < r->seq.end; i++) { _putc(fh, '.'); }            \
      _putc(fh, '\n');                                                         \
    }                                                                          \
    else {                                                                     \
      _printf(fh, "@%s\n", r->name.b);                                         \
      _seq_print_wrap(fh, r->seq.b, r->seq.end, linewrap, i, j, _putc);        \
      _printf(fh, "\n+\n");                                                    \
      _seq_print_wrap(fh, r->qual.b, qlimit, linewrap, i, j, _putc);           \
      /* If i < seq.end, pad quality scores */                                 \
      for(; i < r->seq.end; i++, j++) {                                        \
        if(j == linewrap) { _putc(fh, '\n'); j = 0; }                          \
        _putc(fh, '.');                                                        \
      }                                                                        \
      _putc(fh, '\n');                                                         \
    }                                                                          \
  }

_seq_print_fastq(seq_print_fastq,FILE*,fprintf,fputc2)
_seq_print_fastq(seq_gzprint_fastq,gzFile,gzprintf,gzputc2)

#define SETUP_SEQ_FILE()

// read_t* seq_read_alloc()
// seq_read_destroy(read_t* r)

// seq_open(path)
// seq_open2(path,sam_bam,use_gzip,buffer_size)
// seq_open_fh(fh,use_gzip,buffer_size)
// seq_close(seq_file_t *sf)

#endif
