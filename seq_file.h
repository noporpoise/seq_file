
#ifndef _SEQ_FILE_HEADER
#define _SEQ_FILE_HEADER

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h> // strcasecmp
#include <ctype.h>
#include <limits.h>
#include <zlib.h>

// #define _USESAM

#ifdef _USESAM
  #include "hts.h"
  #include "sam.h"
#endif

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
  #ifdef _USESAM
    samFile *s_file;
    bam_hdr_t *bam_header;
  #endif
  int (*read)(seq_file_t *sf, read_t *r);
  buffer_t in;
  seqtype_t format;
  // Reads pushed onto a 'read stack' aka buffer
  read_t *rhead, *rtail; // 'unread' reads, add to tail, return from head
  int (*origread)(seq_file_t *sf, read_t *r); // used when read = _seq_read_pop
};

struct read_t
{
  buffer_t name, seq, qual;
  #ifdef _USESAM
    bam1_t *bam;
  #endif
  char from_sam; // from sam or bam
  read_t *next; // for use in a linked list
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
#define seq_read(sf,r) ((sf)->read(sf,r))

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

#define DEFAULT_BUFSIZE (1<<20)

//
// Create and destroy read structs
//

static inline void seq_read_dealloc(read_t *r)
{
  if(r->name.b != NULL) free(r->name.b);
  if(r->seq.b != NULL) free(r->seq.b);
  if(r->qual.b != NULL) free(r->qual.b);
  #ifdef _USESAM
    if(r->bam != NULL) free(r->bam);
  #endif
}

static inline read_t* seq_read_alloc(read_t *r)
{
  if(!buffer_init(&r->name, 512) ||
     !buffer_init(&r->seq, DEFAULT_BUFSIZE) ||
     !buffer_init(&r->qual, DEFAULT_BUFSIZE))
  {
    seq_read_dealloc(r);
    return NULL;
  }
  #ifdef _USESAM
    if((r->bam = bam_init1()) == NULL) {
      seq_read_dealloc(r);
      return NULL;
    }
  #endif
  r->from_sam = 0;
  r->next = NULL;
  return r;
}

static inline read_t* seq_read_new()
{
  read_t *r = calloc(1, sizeof(read_t));
  if(r == NULL) return NULL;
  if(seq_read_alloc(r) == NULL) { free(r); return NULL; }
  return r;
}

static inline void seq_read_free(read_t *r)
{
  seq_read_dealloc(r);
  free(r);
}

// file could be sam,bam,FASTA,FASTQ,txt (+gzip)

// Complement SAM/BAM bases
static const int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 9, 14,
                                           1, 6, 5, 13, 3, 11, 7, 15 };


#ifdef _USESAM
// Read a sam/bam file
static inline int _seq_read_sam(seq_file_t *sf, read_t *read)
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
#endif /* _USESAM */


#define _reset_read(sf,read) ({                                                \
  read->name.end = read->seq.end = read->qual.end = 0;                         \
  read->name.b[0] = read->seq.b[0] = read->qual.b[0] = '\0';                   \
  read->from_sam = 0;                                                          \
})


#define _func_read_fastq(_read_fastq,__getc,__ungetc,__readline)               \
  static inline int _read_fastq(seq_file_t *sf, read_t *read)                  \
  {                                                                            \
    int c = __getc(sf);                                                        \
    _reset_read(sf,read);                                                      \
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
    __ungetc(sf, c);                                                           \
    return 1;                                                                  \
  }

#define _func_read_fasta(_read_fasta,__getc,__ungetc,__readline)               \
  static inline int _read_fasta(seq_file_t *sf, read_t *read)                  \
  {                                                                            \
    int c = __getc(sf);                                                        \
    _reset_read(sf,read);                                                      \
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
    __ungetc(sf, c);                                                           \
    return 1;                                                                  \
  }


#define _func_read_plain(_read_plain,__getc,__readline,__skipline)             \
  static inline int _read_plain(seq_file_t *sf, read_t *read)                  \
  {                                                                            \
    int c;                                                                     \
    _reset_read(sf,read);                                                      \
    while((c = __getc(sf)) != -1 && isspace(c)) if(c != '\n') __skipline(sf);  \
    if(c == -1) return 0;                                                      \
    buffer_append_char(&read->seq, c);                                         \
    __readline(sf, read->seq);                                                 \
    buffer_chomp(&(read->seq));                                                \
    return 1;                                                                  \
  }

#define _func_read_unknown(_read_unknown,__getc,__ungetc,__skipline,__fastq,__fasta,__plain)\
  static inline int _read_unknown(seq_file_t *sf, read_t *read)                \
  {                                                                            \
    int c;                                                                     \
    _reset_read(sf, read);                                                     \
    while((c = __getc(sf)) != -1 && isspace(c)) if(c != '\n') __skipline(sf);  \
    if(c == -1) return 0;                                                      \
    if(c == '@') { sf->format = IS_FASTQ; sf->origread = __fastq; }            \
    else if(c == '>') { sf->format = IS_FASTA; sf->origread = __fasta; }       \
    else { sf->format = IS_PLAIN; sf->origread = __plain; }                    \
    __ungetc(sf, c);                                                           \
    return sf->origread(sf,read);                                              \
  }

#define _SF_SWAP(x,y,tmp) ({(tmp) = (x); (x) = (y); (y) = (tmp);})

// Undefined behaviour if you have not previously called _seq_read_shift
static inline int _seq_read_pop(seq_file_t *sf, read_t *read)
{
  buffer_t tmp;
  read_t *next = sf->rhead;
  sf->rhead = sf->rhead->next;
  _SF_SWAP(read->name, next->name, tmp);
  _SF_SWAP(read->seq, next->seq, tmp);
  _SF_SWAP(read->qual, next->qual, tmp);
  seq_read_free(next);
  if(sf->rhead == NULL) {
    sf->read = sf->origread;
    sf->rtail = NULL;
  }
  return 1;
}

// Add a read onto the read buffer (linked list FIFO)
static inline void _seq_read_shift(seq_file_t *sf, read_t *read)
{
  if(sf->rhead == NULL) {
    sf->read = _seq_read_pop;
    sf->rhead = sf->rtail = read;
  }
  else {
    sf->rtail->next = read;
    sf->rtail = read;
  }
  read->next = NULL;
}

// Load reads until we have at least nbases loaded or we hit EOF
static inline void _seq_buffer_reads(seq_file_t *sf, size_t nbases)
{
  int (*read)(seq_file_t *sf, read_t *r) = sf->origread;

  // Sum bases already in buffer
  read_t *r = sf->rhead;
  size_t currbases = 0;
  while(r != NULL) { currbases += r->seq.end; r = r->next; }

  while(currbases < nbases) {
    if((r = seq_read_new()) == NULL) {
      fprintf(stderr, "[%s:%i] Error out of memory\n", __FILE__, __LINE__);
      break;
    }
    if(read(sf,r) <= 0) { seq_read_free(r); break; }
    currbases += r->seq.end;
    _seq_read_shift(sf, r);
  }
}

// perform reading on seq_file_t

// getc on seq_file_t
#define _sf_gzgetc(sf)              gzgetc(sf->gz_file)
#define _sf_gzgetc_buf(sf)          gzgetc_buf(sf->gz_file,&sf->in)
#define _sf_fgetc(sf)               fgetc(sf->f_file)
#define _sf_fgetc_buf(sf)           fgetc_buf(sf->f_file,&sf->in)

// ungetc on seq_file_t
#define _sf_gzungetc(sf,c)          gzungetc(c,sf->gz_file)
#define _sf_gzungetc_buf(sf,c)      ungetc_buf(c,&sf->in)
#define _sf_fungetc(sf,c)           fungetc(c,sf->f_file)
#define _sf_fungetc_buf(sf,c)       ungetc_buf(c,&sf->in)

// readline on seq_file_t using buffer into read
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
_func_read_fastq(_seq_read_fastq_f,      _sf_fgetc,      _sf_fungetc,      _sf_freadline)
_func_read_fastq(_seq_read_fastq_gz,     _sf_gzgetc,     _sf_gzungetc,     _sf_gzreadline)
_func_read_fastq(_seq_read_fastq_f_buf,  _sf_fgetc_buf,  _sf_fungetc_buf,  _sf_freadline_buf)
_func_read_fastq(_seq_read_fastq_gz_buf, _sf_gzgetc_buf, _sf_gzungetc_buf, _sf_gzreadline_buf)

// Read FASTA
_func_read_fasta(_seq_read_fasta_f,      _sf_fgetc,      _sf_fungetc,      _sf_freadline)
_func_read_fasta(_seq_read_fasta_gz,     _sf_gzgetc,     _sf_gzungetc,     _sf_gzreadline)
_func_read_fasta(_seq_read_fasta_f_buf,  _sf_fgetc_buf,  _sf_fungetc_buf,  _sf_freadline_buf)
_func_read_fasta(_seq_read_fasta_gz_buf, _sf_gzgetc_buf, _sf_gzungetc_buf, _sf_gzreadline_buf)

// Read plain
_func_read_plain(_seq_read_plain_f,      _sf_fgetc,      _sf_freadline,      _sf_fskipline)
_func_read_plain(_seq_read_plain_gz,     _sf_gzgetc,     _sf_gzreadline,     _sf_gzskipline)
_func_read_plain(_seq_read_plain_f_buf,  _sf_fgetc_buf,  _sf_freadline_buf,  _sf_fskipline_buf)
_func_read_plain(_seq_read_plain_gz_buf, _sf_gzgetc_buf, _sf_gzreadline_buf, _sf_gzskipline_buf)

// Read first entry
_func_read_unknown(_seq_read_unknown_f,      _sf_fgetc,      _sf_fungetc,      _sf_fskipline,  _seq_read_fastq_f,      _seq_read_fasta_f,      _seq_read_plain_f)
_func_read_unknown(_seq_read_unknown_gz,     _sf_gzgetc,     _sf_gzungetc,     _sf_gzskipline, _seq_read_fastq_gz,     _seq_read_fasta_gz,     _seq_read_plain_gz)
_func_read_unknown(_seq_read_unknown_f_buf,  _sf_fgetc_buf,  _sf_fungetc_buf,  _sf_fskipline,  _seq_read_fastq_f_buf,  _seq_read_fasta_f_buf,  _seq_read_plain_f_buf)
_func_read_unknown(_seq_read_unknown_gz_buf, _sf_gzgetc_buf, _sf_gzungetc_buf, _sf_gzskipline, _seq_read_fastq_gz_buf, _seq_read_fasta_gz_buf, _seq_read_plain_gz_buf)


static inline void _seq_file_init(seq_file_t *sf)
{
  sf->gz_file = NULL;
  sf->f_file = NULL;
  #ifdef _USESAM
    sf->s_file = NULL;
  #endif
  sf->in.size = sf->in.begin = sf->in.end = 0;
  sf->path = sf->in.b = NULL;
  sf->format = IS_UNKNOWN;
  sf->rhead = sf->rtail = NULL;
}

// Returns 1 on success 0 if out of memory
static inline char _seq_setup(seq_file_t *sf, char use_zlib, size_t buf_size)
{
  if(buf_size) {
    if(!buffer_init(&sf->in, buf_size)) { free(sf); return 0; }
    sf->origread = use_zlib ? _seq_read_unknown_gz_buf : _seq_read_unknown_f_buf;
  }
  else sf->origread = use_zlib ? _seq_read_unknown_gz : _seq_read_unknown_f;
  sf->read = sf->origread;
  return 1;
}

// Guess file type from file path or contents
static inline seqtype_t seq_guess_filetype_from_extension(const char *path)
{
  #define _NUM_SEQ_EXT 28
  size_t plen = strlen(path);
  const char *exts[_NUM_SEQ_EXT]
    = {".fa", ".fasta", ".fsa", ".fsa.gz", "fsa.gzip", // FASTA
       ".faz", ".fagz", ".fa.gz", ".fa.gzip", ".fastaz", ".fasta.gzip",
       ".fq", ".fastq", ".fsq", ".fsq.gz", "fsq.gzip", // FASTQ
       ".fqz", ".fqgz", ".fq.gz", ".fq.gzip", ".fastqz", ".fastq.gzip",
       ".txt", ".txtgz", ".txt.gz", ".txt.gzip", // Plain
       ".sam", ".bam"}; // SAM / BAM

  const seqtype_t types[_NUM_SEQ_EXT]
    = {IS_FASTA, IS_FASTA, IS_FASTA, IS_FASTA, IS_FASTA, IS_FASTA,
       IS_FASTA, IS_FASTA, IS_FASTA, IS_FASTA, IS_FASTA,
       IS_FASTQ, IS_FASTQ, IS_FASTQ, IS_FASTQ, IS_FASTQ, IS_FASTQ,
       IS_FASTQ, IS_FASTQ, IS_FASTQ, IS_FASTQ, IS_FASTQ,
       IS_PLAIN, IS_PLAIN, IS_PLAIN, IS_PLAIN,
       IS_SAM, IS_BAM};

  size_t extlens[_NUM_SEQ_EXT];
  size_t i;
  for(i = 0; i < _NUM_SEQ_EXT; i++)
    extlens[i] = strlen(exts[i]);

  for(i = 0; i < _NUM_SEQ_EXT; i++)
    if(extlens[i] <= plen && strcasecmp(path+plen-extlens[i], exts[i]) == 0)
      return types[i];

  return IS_UNKNOWN;
}

// sam_bam is 0 for not SAM/BAM, 1 for SAM, 2 for BAM
static inline seq_file_t* seq_open2(const char *p, char sam_bam,
                                    char use_zlib, size_t buf_size)
{
  seq_file_t *sf = calloc(1, sizeof(seq_file_t));
  _seq_file_init(sf);

  if(sam_bam)
  {
    if(sam_bam != 1 && sam_bam != 2) {
      fprintf(stderr, "[%s:%i] Error: sum_bam param must be 0, 1 (sam) or 2 (bam)\n",
              __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    #ifdef _USESAM
      if((sf->s_file = sam_open(p, sam_bam == 1 ? "rs" : "rb", 0)) == NULL) {
        free(sf->path);
        free(sf);
        return NULL;
      }
      sf->bam_header = sam_hdr_read(sf->s_file);
      sf->read = sf->origread = _seq_read_sam;
      sf->format = sam_bam == 1 ? IS_SAM : IS_BAM;
    #else
      fprintf(stderr, "[%s:%i] Error: not compiled with sam/bam support\n",
              __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    #endif
  }
  else
  {
    if(( use_zlib && ((sf->gz_file = gzopen(p, "r")) == NULL)) ||
       (!use_zlib && ((sf->f_file  =  fopen(p, "r")) == NULL))) {
      free(sf);
      return NULL;
    }

    if(!_seq_setup(sf, use_zlib, buf_size)) return NULL;
  }

  sf->path = strdup(p);
  return sf;
}

// Returns pointer to new seq_file_t on success, seq_close will close the fh,
// so you shouldn't call fclose(fh)
// Returns NULL on error, in which case FILE will not have been closed (caller
// should then call fclose(fh))
static inline seq_file_t* seq_open_fh(FILE *fh, char sam_bam,
                                      char use_zlib, size_t buf_size)
{
  seq_file_t *sf = calloc(1, sizeof(seq_file_t));
  _seq_file_init(sf);

  if(sam_bam)
  {
    if(sam_bam != 1 && sam_bam != 2) {
      fprintf(stderr, "[%s:%i] Error: sum_bam param must be 0, 1 (sam) or 2 (bam)\n",
              __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    #ifdef _USESAM
      if((sf->s_file = sam_open("-", sam_bam == 1 ? "rs" : "rb", 0)) == NULL) {
        free(sf);
        return NULL;
      }
      sf->bam_header = sam_hdr_read(sf->s_file);
      sf->read = sf->origread = _seq_read_sam;
    #else
      fprintf(stderr, "[%s:%i] Error: not compiled with sam/bam support\n",
              __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    #endif
  }
  else
  {
    if(!use_zlib) sf->f_file = fh;
    else if((sf->gz_file = gzdopen(fileno(fh), "r")) == NULL) {
      free(sf);
      return NULL;
    }

    if(!_seq_setup(sf, use_zlib, buf_size)) return NULL;
  }

  sf->path = strdup("-");
  return sf;
}

static inline seq_file_t* seq_open(const char *p)
{
  if(strcmp(p,"-") == 0) return seq_open_fh(stdin, 0, 1, DEFAULT_BUFSIZE);

  seqtype_t format = seq_guess_filetype_from_extension(p);
  char sam_bam = 0;
  if(format == IS_SAM) sam_bam = 1;
  if(format == IS_BAM) sam_bam = 2;
  return seq_open2(p, sam_bam, 1, DEFAULT_BUFSIZE);
}

static inline void seq_close(seq_file_t *sf)
{
  if(sf->f_file != NULL) { fclose(sf->f_file); }
  if(sf->gz_file != NULL) { gzclose(sf->gz_file); }
  #ifdef _USESAM
  if(sf->s_file != NULL) { sam_close(sf->s_file); free(sf->bam_header); }
  #endif
  if(sf->in.b != NULL) { free(sf->in.b); }
  if(sf->path != NULL) { free(sf->path); }
  read_t *r = sf->rhead, *tmpr;
  while(r != NULL) { tmpr = r->next; seq_read_free(r); r = tmpr; }
  free(sf);
}

// Get min and max quality values by examining the first 
// Returns 0 if no qual scores, 1 on success, -1 if read error
static inline int seq_get_qual_limits(seq_file_t *sf, int *minq, int *maxq)
{
  read_t *r;
  int min = INT_MAX, max = 0;
  size_t i, count = 0, qcount = 0, limit = 1000, len;
  char q;

  _seq_buffer_reads(sf, limit);
  r = sf->rhead;

  while(count < limit && r != NULL)
  {
    len = (r->qual.end < limit - qcount ? r->qual.end : limit - qcount);
    for(i = 0; i < len; i++) {
      q = r->qual.b[i];
      if(q > max) max = q;
      if(q < min) min = q;
    }
    count += r->seq.end;
    qcount += r->qual.end;
    r = r->next;
  }

  if(qcount > 0) { *minq = min; *maxq = max; }

  return (qcount > 0);
}

// Returns -1 on error
// Returns 0 if not in FASTQ format/ not recognisable (offset:33, min:33, max:104)
// max_read_bases is the number of quality scores to read
static inline int seq_guess_fastq_format(seq_file_t *sf, int *minq, int *maxq)
{
  // Detect fastq offset
  *minq = INT_MAX;
  *maxq = 0;

  if(seq_get_qual_limits(sf, minq, maxq) <= 0) return -1;

  // See: http://en.wikipedia.org/wiki/FASTQ_format
  // Usually expect 0,40, but new software can report 41, so using <= MAX+1
  if(*minq >= 33 && *maxq <= 73) return 1; // sanger
  else if(*minq >= 33 && *maxq <= 75) return 5; // Illumina 1.8+
  else if(*minq >= 67 && *maxq <= 105) return 4; // Illumina 1.5+
  else if(*minq >= 64 && *maxq <= 105) return 3; // Illumina 1.3+
  else if(*minq >= 59 && *maxq <= 105) return 2; // Solexa
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
#define seq_read_looks_valid_protein(r) \
        _seq_read_looks_valid(r,"acdefghiklmnopqrstuvwy")

static inline char _seq_char_complement(char c) {
  switch(c) {
    case 'a': return 't'; case 'A': return 'T';
    case 'c': return 'g'; case 'C': return 'G';
    case 'g': return 'c'; case 'G': return 'C';
    case 't': return 'a'; case 'T': return 'A';
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

// Formally, FASTA/Q entry names stop at the first space character
// Truncates read name and returns new length
static inline size_t seq_read_truncate_name(read_t *r)
{
  char *tmp = r->name.b;
  size_t len;
  for(len = 0; tmp[len] != '\0' && !isspace(tmp[len]); len++) {}
  r->name.b[r->name.end = len] = '\0';
  return len;
}

static inline void seq_read_to_uppercase(read_t *r)
{
  char *tmp;
  for(tmp = r->seq.b; *tmp != '\0'; tmp++) *tmp = toupper(*tmp);
}

static inline void seq_read_to_lowercase(read_t *r)
{
  char *tmp;
  for(tmp = r->seq.b; *tmp != '\0'; tmp++) *tmp = tolower(*tmp);
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

// New read on the stack
// read_t* seq_read_alloc(read_t*)
// seq_read_dealloc(read_t* r)

// New read on the heap
// read_t* seq_read_new()
// seq_read_free(read_t* r)

// seq_open(path)
// seq_open2(path,sam_bam,use_gzip,buffer_size)
// seq_open_fh(fh,use_gzip,buffer_size)
// seq_close(seq_file_t *sf)

#endif
